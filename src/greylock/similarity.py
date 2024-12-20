"""Module for calculating weighted subcommunity and metacommunity 
similarities.

Classes
-------
Similarity
    Abstract base class for relative abundance-weighted species
    similarities.
SimilarityFromDataFrame
    Implements Similarity by storing similarities in a pandas DataFrame.
SimilarityFromArray
    Implements Similarity by storing similarities in a numpy ndarray or 
    memmap.
SimilarityFromFile
    Implements Similarity by reading similarities from a csv or tsv 
    file.
SimilarityFromFunction:
    Implements Similarity by calculating pairwise similarities with a 
    callable function.

Functions
---------
weighted_similarity_chunk_nonsymmetric
     Utility funciton for calculation of chunks of a similarity matrix
weighted_similarity_chunk_symmetric
     Utility funciton for calculation of chunks of a symmetric similarity matrix 
"""

from abc import ABC, abstractmethod
from typing import Callable, Union, Tuple
from pathlib import Path
from numpy import ndarray, empty, concatenate, float64, vstack, zeros
from pandas import DataFrame, read_csv
from scipy.sparse import spmatrix, issparse  # type: ignore[import]
from greylock.exceptions import InvalidArgumentError


class Similarity(ABC):
    """Interface for classes computing weighted similarities."""

    @abstractmethod
    def weighted_abundances(
        self, relative_abundances: Union[ndarray, spmatrix]
    ) -> ndarray:
        """Calculates weighted sums of similarities to each species.

        Returns
        -------
        A 2-d numpy.ndarray of shape (n_species, n_communities), where
        rows correspond to unique species, columns correspond to
        (meta-/sub-) communities and each element is a sum of
        similarities of all species to the row's species weighted by
        their relative abundances in the respective communities.
        """
        pass

    def self_similar_weighted_abundances(
        self, relative_abundances: Union[ndarray, spmatrix]
    ) -> ndarray:
        return self.weighted_abundances(relative_abundances)

    def is_expensive(self):
        return False

    def __matmul__(self, abundance):
        return abundance.premultiply_by(self)


class SimilarityIdentity(Similarity):
    def weighted_abundances(self, relative_abundance):
        return relative_abundance


class SimilarityFromArray(Similarity):
    """Implements Similarity using similarities stored in a numpy
    ndarray."""

    def __init__(self, similarity: Union[ndarray, spmatrix]):
        self.similarity = similarity

    def weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ) -> ndarray:
        return self.similarity @ relative_abundance

    def self_similar_weighted_abundances(
        self, relative_abundances: Union[ndarray, spmatrix]
    ) -> ndarray:
        if self.similarity.shape[0] == self.similarity.shape[1]:
            return self.weighted_abundances(relative_abundances)
        else:
            raise InvalidArgumentError("Similarity matrix must be square")


class SimilarityFromDataFrame(SimilarityFromArray):
    """Implements Similarity using similarities stored in pandas
    dataframe."""

    def __init__(self, similarity: DataFrame):
        self.similarity = similarity

    def weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ) -> ndarray:
        return self.similarity.to_numpy() @ relative_abundance


class SimilarityFromFile(Similarity):
    """Implements Similarity by using similarities stored in file.

    Similarity matrix rows are read from the file one chunk at a time.
    The size of chunks can be specified in numbers of rows to control
    memory load.
    """

    def __init__(
        self, similarity_file_path: Union[str, Path], chunk_size: int = 100
    ) -> None:
        """
        Parameters
        ----------
        similarity_file_path:
            Path to a file containing a pairwise similarity matrix of
            shape (n_species, n_species). The file should have a header
            that denotes the unique species names.
        chunk_size:
            Number of rows to read from similarity matrix at a time.
        """
        self.path = Path(similarity_file_path)
        self.chunk_size = chunk_size

    def weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ) -> ndarray:
        weighted_abundances = empty(relative_abundance.shape, dtype=float64)
        with read_csv(
            self.path,
            chunksize=self.chunk_size,
            sep=None,
            engine="python",
            dtype=float64,
        ) as similarity_matrix_chunks:
            i = 0
            for chunk in similarity_matrix_chunks:
                weighted_abundances[i : i + self.chunk_size, :] = (
                    chunk.to_numpy() @ relative_abundance
                )
                i += self.chunk_size
        return weighted_abundances

    def is_expensive(self):
        return True


class IntersetSimilarityFromFile(SimilarityFromFile):
    def weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ) -> ndarray:
        with read_csv(
            self.path,
            chunksize=self.chunk_size,
            sep=None,
            engine="python",
            dtype=float64,
        ) as similarity_matrix_chunks:
            weighted_abundance_chunks = [
                chunk.to_numpy() @ relative_abundance
                for chunk in similarity_matrix_chunks
            ]
        return concatenate(weighted_abundance_chunks)

    def self_similar_weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ):
        raise InvalidArgumentError(
            "Inappropriate similarity class for diversity measures"
        )


def weighted_similarity_chunk_nonsymmetric(
    similarity: Callable,
    X: Union[ndarray, DataFrame],
    Y: Union[ndarray, DataFrame, None],
    relative_abundance: ndarray,
    chunk_size: int,
    chunk_index: int,
) -> Tuple[int, ndarray]:
    def enum_helper(X):
        if type(X) == DataFrame:
            return X.itertuples()
        return X

    if Y is None:
        Y = X
    chunk = X[chunk_index : chunk_index + chunk_size]
    similarities_chunk = empty(shape=(chunk.shape[0], Y.shape[0]))
    for i, row_i in enumerate(enum_helper(chunk)):
        for j, row_j in enumerate(enum_helper(Y)):
            similarities_chunk[i, j] = similarity(row_i, row_j)
    # When this is a remote task, the chunks may be returned out of
    # order. Indicate what chunk this was for, so we can sort the
    # resulting chunks correctly:
    return chunk_index, similarities_chunk @ relative_abundance


def weighted_similarity_chunk_symmetric(
    similarity: Callable,
    X: Union[ndarray, DataFrame],
    relative_abundance: ndarray,
    chunk_size: int,
    chunk_index: int,
) -> ndarray:
    def enum_helper(X, start_index=0):
        if type(X) == DataFrame:
            return X.iloc[start_index:].itertuples()
        return X[start_index:]

    chunk = X[chunk_index : chunk_index + chunk_size]
    similarities_chunk = zeros(shape=(chunk.shape[0], X.shape[0]))
    for i, row_i in enumerate(enum_helper(chunk)):
        for j, row_j in enumerate(enum_helper(X, chunk_index + i + 1)):
            similarities_chunk[i, i + j + chunk_index + 1] = similarity(row_i, row_j)
    rows_result = similarities_chunk @ relative_abundance
    rows_after_count = max(0, relative_abundance.shape[0] - (chunk_index + chunk_size))
    rows_result = vstack(
        (
            zeros(shape=(chunk_index, relative_abundance.shape[1])),
            rows_result,
            zeros(
                shape=(
                    rows_after_count,
                    relative_abundance.shape[1],
                )
            ),
        )
    )
    similarities_chunk = similarities_chunk.T
    relative_abundance = relative_abundance[chunk_index : chunk_index + chunk_size]
    cols_result = similarities_chunk @ relative_abundance
    return rows_result + cols_result


class SimilarityFromSymmetricFunction(Similarity):
    """
    Calculate a similarity matrix on the fly, given feature vectors for each
    species and a function to calculate similarity from feature vectors.

    This assumes that
    a) the similarity between any species and itself is 1.0.
    b) similarity is symmetric: that is, Z[i,j] == Z[j,i]] always by definition.

    Note that condition (b) need not be true and all the Leinster & Cobbald math
    is still valid (see L&C paper). If these conditions do not hold, use the
    next class: SimilarityFromFunction. But if you can use this class, you get the
    obvious speed-up.

    N.B. that these calculations can be parallelized; see ray.py for
    parallelization using the Ray package.
    """

    def __init__(
        self, func: Callable, X: Union[ndarray, DataFrame], chunk_size: int = 100
    ):
        """
        Parameters
        ----------
        func:
            A Callable that calculates similarity between a pair of species.
            Must take two rows from X as input as its arguments, and return
            a numeric similarity value between 0.0 and 1.0.
        X:
          Each row contains the features values for a given species.
        chunk_size:
            Number of rows in similarity matrix to calculate at a time.
        """
        self.func = func
        self.X = X
        self.chunk_size = chunk_size

    def is_expensive(self):
        return True

    def weighted_abundances(self, abundance: Union[ndarray, spmatrix]) -> ndarray:
        result = abundance.copy()
        for chunk_index in range(0, self.X.shape[0], self.chunk_size):
            chunk = weighted_similarity_chunk_symmetric(
                similarity=self.func,
                X=self.X,
                relative_abundance=abundance,
                chunk_size=self.chunk_size,
                chunk_index=chunk_index,
            )
            result = result + chunk
        return result


class SimilarityFromFunction(SimilarityFromSymmetricFunction):
    def weighted_abundances(self, relative_abundance: Union[ndarray, spmatrix]):
        weighted_similarity_chunks = []
        for chunk_index in range(0, self.X.shape[0], self.chunk_size):
            _, result = weighted_similarity_chunk_nonsymmetric(
                self.func,
                self.X,
                self.get_Y(),
                relative_abundance,
                self.chunk_size,
                chunk_index,
            )
            weighted_similarity_chunks.append(result)
        return concatenate(weighted_similarity_chunks)

    def get_Y(self):
        return None


class IntersetSimilarityFromFunction(SimilarityFromFunction):
    def __init__(
        self,
        func: Callable,
        X: Union[ndarray, DataFrame],
        Y: Union[ndarray, DataFrame],
        chunk_size: int = 100,
    ):
        """
        Parameters
        ----------
        func:
            A Callable that calculates similarity between a pair of species.
            Must take a row from X and a row from Y as input as its arguments, and return
            a numeric similarity value between 0.0 and 1.0.
        X:
          Each row contains the features values for a given species in set A.
        Y:
          Each row contains the features values for a given species in set B.
        chunk_size:
            Number of rows in similarity matrix to calculate at a time.
        """
        super().__init__(func, X, chunk_size)
        self.Y = Y

    def get_Y(self):
        return self.Y

    def self_similar_weighted_abundances(
        self, relative_abundance: Union[ndarray, spmatrix]
    ):
        raise InvalidArgumentError(
            "Inappropriate similarity class for diversity measures"
        )
