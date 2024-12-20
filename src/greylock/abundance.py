"""Module for calculating relative sub- and metacomunity abundances.

Classes
-------
Abundance
    Relative (normalized) species abundances in (meta-/sub-) communities
AbundanceForDiversity
    Species abundances-- normalized over metacommunity, normalized over each subcommunity,
    and totalled across metacommunity-- as is required for diversity calculations

"""

from functools import cached_property
from typing import Iterable, Union

from numpy import arange, ndarray, concatenate
from pandas import DataFrame, RangeIndex
from scipy.sparse import issparse  # type: ignore[import]


class Abundance:
    def __init__(
        self, counts: ndarray, subcommunity_names: Iterable[Union[str, int]]
    ) -> None:
        """
        Parameters
        ----------
        counts:
            2-d array with one column per subcommunity, one row per
            species, containing the count of each species in the
            corresponding subcommunities.
        """
        self.subcommunities_names = subcommunity_names
        self.num_subcommunities = counts.shape[1]
        self.subcommunity_abundance = self.make_subcommunity_abundance(counts=counts)
        self.normalized_subcommunity_abundance = (
            self.make_normalized_subcommunity_abundance()
        )

    def make_subcommunity_abundance(self, counts: ndarray) -> ndarray:
        """Calculates the relative abundances in subcommunities.

        Parameters
        ----------
        counts
            2-d array with one column per subcommunity, one row per
            species, containing the count of each species in the
            corresponding subcommunities.

        Returns
        -------
        A numpy.ndarray of shape (n_species, n_subcommunities), where
        rows correspond to unique species, columns correspond to
        subcommunities and each element is the abundance of the species
        in the subcommunity relative to the total metacommunity size.
        """
        return counts / counts.sum()

    def make_subcommunity_normalizing_constants(self) -> ndarray:
        """Calculates subcommunity normalizing constants.

        Returns
        -------
        A numpy.ndarray of shape (n_subcommunities,), with the fraction
        of each subcommunity's size of the metacommunity.
        """
        return self.subcommunity_abundance.sum(axis=0)

    def make_normalized_subcommunity_abundance(self) -> ndarray:
        """Calculates normalized relative abundances in subcommunities.

        Returns
        -------
        A numpy.ndarray of shape (n_species, n_subcommunities), where
        rows correspond to unique species, columns correspond to
        subcommunities and each element is the abundance of the species
        in the subcommunity relative to the subcommunity size.
        """
        self.subcommunity_normalizing_constants = (
            self.make_subcommunity_normalizing_constants()
        )
        return self.subcommunity_abundance / self.subcommunity_normalizing_constants

    def premultiply_by(self, similarity):
        return similarity.weighted_abundances(self.normalized_subcommunity_abundance)


class AbundanceForDiversity(Abundance):
    """Calculates metacommuntiy and subcommunity relative abundance
    components from a numpy.ndarray containing species counts
    """

    def __init__(
        self, counts: ndarray, subcommunity_names: Iterable[Union[str, int]]
    ) -> None:
        super().__init__(counts, subcommunity_names)
        self.metacommunity_abundance = self.make_metacommunity_abundance()
        self.unified_abundance_array = None

    def unify_abundance_array(self) -> None:
        """Creates one matrix containing all the abundance matrices:
        metacommunity, subcommunity, and normalized subcommunity.
        These matrices are still available as views on the unified
        data structure. (Because we are using basic slicing here, only
        one copy of the data will exist after garbage collection.)

        This allows for a major computational improvement in efficiency:
        The similarity matrix only has to be generated and used
        once (in the case where a pre-computed similarity matrix is not
        in RAM). That is, we make only one call to
        similarity.weighted_abundances(), in cases where generation of the
        similarity matrix is expensive.
        """
        self.unified_abundance_array = concatenate(
            (
                self.metacommunity_abundance,
                self.subcommunity_abundance,
                self.normalized_subcommunity_abundance,
            ),
            axis=1,
        )

    def get_unified_abundance_array(self):
        if self.unified_abundance_array is None:
            self.unify_abundance_array()
            self.metacommunity_abundance = self.unified_abundance_array[:, [0]]
            self.subcommunity_abundance = self.unified_abundance_array[
                :, 1 : (1 + self.num_subcommunities)
            ]
            self.normalized_subcommunity_abundance = self.unified_abundance_array[
                :, (1 + self.num_subcommunities) :
            ]
        return self.unified_abundance_array

    def premultiply_by(self, similarity):
        if similarity.is_expensive():
            all_ordinariness = similarity.self_similar_weighted_abundances(
                self.get_unified_abundance_array()
            )
            metacommunity_ordinariness = all_ordinariness[:, [0]]
            subcommunity_ordinariness = all_ordinariness[
                :, 1 : (1 + self.num_subcommunities)
            ]
            normalized_subcommunity_ordinariness = all_ordinariness[
                :, (1 + self.num_subcommunities) :
            ]
        else:
            metacommunity_ordinariness = similarity.self_similar_weighted_abundances(
                self.metacommunity_abundance
            )
            subcommunity_ordinariness = similarity.self_similar_weighted_abundances(
                self.subcommunity_abundance
            )
            normalized_subcommunity_ordinariness = (
                similarity.self_similar_weighted_abundances(
                    self.normalized_subcommunity_abundance
                )
            )
        return (
            metacommunity_ordinariness,
            subcommunity_ordinariness,
            normalized_subcommunity_ordinariness,
        )

    def make_metacommunity_abundance(self) -> ndarray:
        """Calculates the relative abundances in metacommunity.

        Returns
        -------
        A numpy.ndarray of shape (n_species, 1), where rows correspond
        to unique species and each row contains the relative abundance
        of the species in the metacommunity.
        """
        return self.subcommunity_abundance.sum(axis=1, keepdims=True)


def make_abundance(counts: Union[DataFrame, ndarray], for_diversity=True) -> Abundance:
    """Initializes a concrete subclass of Abundance.

    Parameters
    ----------
    counts:
        2-d array with one column per subcommunity, one row per species,
        where the elements are the species counts

    Returns
    -------
    An instance of a concrete subclass of Abundance.
    """
    if for_diversity:
        specific_class = AbundanceForDiversity
    else:
        specific_class = Abundance
    if isinstance(counts, DataFrame):
        return specific_class(
            counts=counts.to_numpy(), subcommunity_names=counts.columns
        )
    elif hasattr(counts, "shape"):
        if issparse(counts):
            raise TypeError("sparse abundance matrix not yet implemented")
        else:
            return specific_class(
                counts=counts, subcommunity_names=arange(counts.shape[1])
            )
    else:
        raise NotImplementedError(
            f"Type {type(counts)} is not supported for argument "
            "'counts'. Valid types include pandas.DataFrame or"
            "numpy.ndarray"
        )
