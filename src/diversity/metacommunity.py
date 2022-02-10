"""Module for metacommunity and subcommunity diversity measures.

Classes
-------
Metacommunity
    Represents a metacommunity made up of subcommunities and computes
    metacommunity subcommunity diversity measures.

Functions
---------
make_metacommunity
    Builds diversity.metacommunity.Metacommunity object according to
    parameter specification.
"""

from abc import ABC, abstractmethod
from functools import cache

from pandas import DataFrame, concat, unique
from numpy import zeros, broadcast_to, divide

from diversity.abundance import Abundance
from diversity.log import LOGGER
from diversity.similarity import make_similarity
from diversity.utilities import (
    pivot_table,
    power_mean,
    subset_by_column,
)

def make_metacommunity(
    counts,
    similarity_method,
    subcommunities=None,
    shared_abundance=False,
    chunk_size=1,
    features_filepath=None,
    shared_array_manager=None,
    num_processors=None,
    subcommunity_column="subcommunity",
    species_column="species",
    count_column="count",
):
    """Builds a Metacommunity object from specified parameters.

    Depending on the chosen counts argument, the subcommunity and
    normalized subcommunity relative abundances are either
    simultaneously stored in memory for memory-heavy but fast
    computations (when using a pandas.DataFrame), or stored one at a
    time in a shared memory block for large data sets (when providing a
    filepath).

    The similarity_method parameter value determines whether:
        1. An in-memory similarity matrix is used, or
        2. Similarities are read from a file, or
        3. Similarities are computed on the fly without storing the
           entire species similarity matrix in memory or on disk.

    Parameters
    ----------
    counts: pandas.DataFrame, or str
        Table or path to file containing table with 3 columns: one
        column lists subcommunity identifiers, one lists species
        identifiers, and the last lists counts of species in
        corresponding subcommunities. Subcommunity-species identifier
        pairs are assumed to be unique. If using Callable as
        similarity_method, this must be a filepath. Both file and data
        frame are assumed to have column headers as specified by the
        subcommunity_column, species_column, and count_column arguments.
    similarity_method: pandas.DataFrame, str, or Callable
        For in-memory data frame, see diversity.similarity.SimilarityFromMemory,
        for str filepath to similarity matrix, see
        diversity.similarity.SimilarityFromFile, and for callable which
        calculates similarities on the fly,
        see diversity.similarity.SimilarityFromFunction.
    subcommunities: collection of str objects supporting membership test
        Names of subcommunities to include. Their union is the
        metacommunity, and data for all other subcommunities is ignored.
    chunk_size: int
        See diversity.similarity.SimilarityFromFile. Only
        relevant when using in-file similarity matrix (str filepath as
        similarity_method).
    features_filepath: str
        Path to .tsv, or .csv file containing species features. Assumed
        to have a header row, with species identifiers in a column with
        header species_column. All other columns are assumed to be
        features. This parameter is only relevant when similarities are
        computed on the fly (when a Callable is used as the
        similarity_method argument).
    num_processors: int
        See diversity.similarity.SimilarityFromFunction. Only relevant
        when calculating similarities on the fly (when a Callable is
        used as the similarity_method argument).
    subcommunity_column, species_column, count_column: str
        Used to specify non-default column headers in counts table.

    Returns
    -------
    A diversity.metacommunity.Metacommunity object built according to
    parameter specification.
    """
    LOGGER.debug(
        "make_metacommunity(counts=%s, similarity_method=%s,"
        " subcommunities=%s, chunk_size=%s, features_filepath=%s,"
        " num_processors=%s, subcommunity_column=%s, species_column=%s,"
        " count_column=%s",
        counts,
        similarity_method,
        subcommunities,
        chunk_size,
        features_filepath,
        num_processors,
        subcommunity_column,
        species_column,
        count_column,
    )

    counts_subset = subset_by_column(counts, subcommunities, subcommunity_column)
    species_subset = unique(counts_subset[species_column])
    similarity = make_similarity(similarity_matrix, species_subset, chunk_size)
    abundance = Abundance(
        pivot_table(
            data_frame=counts_subset,
            pivot_column=subcommunity_column,
            index_column=species_column,
            value_columns=[count_column],
            index_ordering=similarity.species_order,
        )
    )
    return Metacommunity(similarity, abundance)


def make_pairwise_metacommunities(
    counts, similarity_matrix, subcommunity_column, **kwargs
):
    subcommunties_groups = counts.groupby(subcommunity_column)
    pairwise_metacommunities = []
    for i, (_, group_i) in enumerate(subcommunties_groups):
        for j, (_, group_j) in enumerate(subcommunties_groups):
            if j > i:
                counts = concat([group_i, group_j])
                pair_ij = make_metacommunity(
                    counts,
                    similarity_matrix,
                    subcommunity_column=subcommunity_column,
                    **kwargs,
                )
                pairwise_metacommunities.append(pair_ij)
    return pairwise_metacommunities


class IMetacommunity(ABC):
    """Interface for metacommunities and calculating their diversity."""

    @abstractmethod
    def __init__(self, abundance):
        self.__abundance = abundance
        self.__measure_components = None

    @cache
    def subcommunity_diversity(self, viewpoint, measure):
        """Calculates subcommunity diversity measures.

        Parameters
        ----------
        viewpoint: numeric
            Viewpoint parameter for diversity measure.
        measure: str
            Name of the diversity measure.

        Returns
        -------
        A numpy array with a diversity value per subcommunity.

        Notes
        -----
        Valid measure identifiers are: "alpha", "rho", "beta", "gamma",
        "normalized_alpha", "normalized_rho", and "normalized_beta".
        """
        numerator, denominator = self.__measure_components[measure]
        if callable(numerator):
            numerator = numerator()
        denominator = denominator()
        if measure == "gamma":
            denominator = broadcast_to(
                denominator,
                self.__abundance.normalized_subcommunity_abundance().shape,
            )
        community_ratio = divide(
            numerator, denominator, out=zeros(denominator.shape), where=denominator != 0
        )
        result = power_mean(
            1 - viewpoint,
            self.__abundance.normalized_subcommunity_abundance()
            community_ratio,
        )
        if measure in ["beta", "normalized_beta"]:
            return 1 / result
        return result

    @cache
    def metacommunity_diversity(self, viewpoint, measure):
        """Calculates metcommunity diversity measures."""
        subcommunity_diversity = self.subcommunity_diversity(viewpoint, measure)
        return power_mean(
            1 - viewpoint,
            self.__abundance.subcommunity_normalizing_constants()
            subcommunity_diversity,
        )

    def subcommunities_to_dataframe(self, viewpoint):
        """Table containing all subcommunity diversity values.

        Parameters
        ----------
        viewpoint: numeric
            Non-negative number. Can be interpreted as the degree of
            ignorance towards rare species, where 0 treats rare species
            the same as frequent species, and infinity considers only the
            most frequent species.
        """
        df = DataFrame(
            {
                key: self.subcommunity_diversity(viewpoint, key)
                for key in self.__measure_components.keys()
            }
        )
        df.insert(0, "viewpoint", viewpoint)
        df.insert(0, "community", self.__abundance.subcommunity_order)
        return df

    def metacommunity_to_dataframe(self, viewpoint):
        """Table containing all metacommunity diversity values.

        Parameters
        ----------
        viewpoint: numeric
            Non-negative number. Can be interpreted as the degree of
            ignorance towards rare species, where 0 treats rare species
            the same as frequent species, and infinity considers only the
            most frequent species.
        """
        df = DataFrame(
            {
                key: self.metacommunity_diversity(viewpoint, key)
                for key in self.__measure_components.keys()
            },
            index=["metacommunity"],
        )
        df.insert(0, "viewpoint", viewpoint)
        return df


class ISimilaritySensitiveMetacommunity(IMetacommunity):
    """Interface for calculating similarity-sensitive diversity."""

    def __init__(self, abundance, similarity):
        """Initializes object.

        Parameters
        ----------
        abundance: diversity.abundance.IAbundance
            Object whose (sub-/meta-)community species abundances are
            used.
        similarity: diversity.similarity.ISimilarity
            Object for calculating abundance-weighted similarities.
        """
        super().__init__(abundance)
        self.__similarity = similarity
        self.__measure_components = {
            "alpha": (1, self.subcommunity_similarity),
            "rho": (self.metacommunity_similarity, self.subcommunity_similarity),
            "beta": (self.metacommunity_similarity, self.subcommunity_similarity),
            "gamma": (1, self.metacommunity_similarity),
            "normalized_alpha": (1, self.normalized_subcommunity_similarity),
            "normalized_rho": (
                self.metacommunity_similarity,
                self.normalized_subcommunity_similarity,
            ),
            "normalized_beta": (
                self.metacommunity_similarity,
                self.normalized_subcommunity_similarity,
            ),
        }

    @abstractmethod
    def metacommunity_similarity(self):
        """Sums of similarities weighted by metacommunity abundances."""
        pass

    @abstractmethod
    def subcommunity_similarity(self):
        """Sums of similarities weighted by subcommunity abundances."""
        pass

    @abstractmethod
    def normalized_subcommunity_similarity(self):
        """Sums of similarities weighted by the normalized subcommunity abundances."""
        pass


class SimilarityInsensitiveMetacommunity(IMetacommunity):
    """Implements IMetacommunity for similarity-insensitive diversity."""

    def __init__(self, abundance):
        """Initializes object.

        Parameters
        ----------
        abundance: diversity.abundance.IAbundance
            Object whose (sub-/meta-)community species abundances are
            used.
        """
        super().__init__(abundance)
        self.__measure_components = {
            "alpha": (1, self.__abundance.subcommunity_abundance),
            "rho": (
                self.__abundance.metacommunity_abundance,
                self.__abundance.subcommunity_abundance,
            ),
            "beta": (
                self.__abundance.metacommunity_abundance,
                self.__abundance.subcommunity_abundance,
            ),
            "gamma": (1, self.__abundance.metacommunity_abundance),
            "normalized_alpha": (
                1,
                self.__abundance.normalized_subcommunity_abundance,
            ),
            "normalized_rho": (
                self.__abundance.metacommunity_abundance,
                self.__abundance.normalized_subcommunity_abundance,
            ),
            "normalized_beta": (
                self.__abundance.metacommunity_abundance,
                self.__abundance.normalized_subcommunity_abundance,
            ),
        }


class SimilaritySensitiveMetacommunity(ISimilaritySensitiveMetacommunity):
    """Implements ISimilaritySensitiveMetacommunity for fast but memory heavy calculations."""

    @cache
    def metacommunity_similarity(self):
        return self.__similarity.calculate_weighted_similarities(
            self.__abundance.metacommunity_abundance()
        )

    @cache
    def subcommunity_similarity(self):
        return self.__similarity.calculate_weighted_similarities(
            self.__abundance.subcommunity_abundance()
        )

    @cache
    def normalized_subcommunity_similarity(self):
        return self.__similarity.calculate_weighted_similarities(
            self.__abundance.normalized_subcommunity_abundance()
        )

class SharedSimilaritySensitiveMetacommunity(ISimilaritySensitiveMetacommunity):
    """Implements ISimilaritySensitiveMetacommunity using shared memory.

    Caches only one of weighted subcommunity similarities and normalized
    weighted subcommunity similarities at a time. All weighted similarities
    are stored in shared arrays, which can be passed to other processors
    without copying.
    """

    def __init__(self, abundance, similarity, shared_array_manager):
        """Initializes object.

        Parameters
        ----------
        abundance, similarity
            See diversity.metacommunity.ISimilaritySensitiveMetacommunity.
        shared_memory_manager: diversity.shared.SharedMemoryManager
            Active manager for obtaining shared arrays.

        Notes
        -----
        - Object will break once shared_array_manager becomes inactive.
        - If a diversity.similarity.SimilarityFromFunction object is
          chosen as argument for simlarity, it must be paired with
          diversity.abundance.SharedAbundance object as argument for
          abundance.
        """
        super().__init__(abundance=abundance, similarity=similarity)
        self.__storing_normalized_similarities = None
        self.__shared_array_manager = shared_array_manager
        self.__shared_similarity = self.__shared_array_manager.empty(
            shape=self.__abundance.subcommunity_abundance().shape,
            data_type=self.__abundance.subcommunity_abundance().dtype,
        )
        self.__metacommunity_similarity = None


    def metacommunity_similarity(self):
        self.__metacommunity_similarity is None:
            self.__metacommunity_similarity = self.__shared_array_manager.empty(
                shape=self.__abundance.metacommunity_abundance().shape,
                dtype=self.__abundance.metacommunity_abundance().dtype,
            )
            self.__similarity.calculate_weighted_similarities(
                self.__abundance.metacommunity_abundance(),
                out=self.__metacommunity_similarity
            )
        return self.__metacommunity_similarity.data

    @cache
    def subcommunity_similarity(self):
        if self.__storing_normalized_similarities is None:
            self.__similarity.calculate_weighted_similarities(
                self.__abundance.subcommunity_abundance(),
                out=self.__shared_similarity
            )
        elif self.__storing_normalized_similarities:
            self.__shared_similarity.data *= self.__abundance.subcommunity_normalizing_constants()
        self.__storing_normalized_similarities = False
        return self.__shared_similarity.data

    @cache
    def normalized_subcommunity_similarity(self):
        if self.__storing_normalized_similarities is None:
            self.__similarity.calculate_weighted_similarities(
                self.__abundance.normalized_subcommunity_abundance(),
                out=self.__shared_similarity
            )
        elif not self.__storing_normalized_similarities:
            self.__shared_similarity.data /= self.__abundance.subcommunity_normalizing_constants()
        self.__storing_normalized_similarities = True
        return self.__shared_similarity.data
