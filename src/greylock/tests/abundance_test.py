"""Tests for diversity.abundance."""

from dataclasses import dataclass, field
from numpy import allclose, array, ndarray
from pandas import DataFrame
from pytest import fixture, mark, raises
from scipy.sparse import coo_array

from greylock.abundance import (
    AbundanceForDiversity,
    make_abundance,
)


def test_sparse():
    sparse_abundance = coo_array(
        (array([4, 5, 7, 9]), (array([0, 3, 1, 0]), array([0, 3, 1, 2]))),
        shape=(4, 4),
    )
    with raises(TypeError):
        make_abundance(sparse_abundance)


def counts_array_3by2():
    return array([[2, 4], [3, 0], [0, 1]])


counts_dataframe_3by2 = DataFrame(counts_array_3by2())


@dataclass
class AbundanceExclusiveSpecies:
    description: str = "2 subcommunities; both contain exclusive species"
    counts: ndarray = field(default_factory=counts_array_3by2)
    subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [[2 / 10, 4 / 10], [3 / 10, 0 / 10], [0 / 10, 1 / 10]]
        )
    )
    metacommunity_abundance: ndarray = field(
        default_factory=lambda: array([[6 / 10], [3 / 10], [1 / 10]])
    )
    subcommunity_normalizing_constants: ndarray = field(
        default_factory=lambda: array([5 / 10, 5 / 10])
    )
    normalized_subcommunity_abundance: ndarray = field(
        default_factory=lambda: array([[2 / 5, 4 / 5], [3 / 5, 0 / 5], [0 / 5, 1 / 5]])
    )


@dataclass
class AbundanceSingleExclusiveSpecies:
    description: str = "2 subcommunities; one contains exclusive species"
    counts: ndarray = field(default_factory=lambda: array([[2, 4], [3, 0], [5, 1]]))
    subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [[2 / 15, 4 / 15], [3 / 15, 0 / 15], [5 / 15, 1 / 15]]
        )
    )
    metacommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [6 / 15],
                [3 / 15],
                [6 / 15],
            ]
        )
    )
    subcommunity_normalizing_constants: ndarray = field(
        default_factory=lambda: array([10 / 15, 5 / 15])
    )
    normalized_subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [[2 / 10, 4 / 5], [3 / 10, 0 / 5], [5 / 10, 1 / 5]]
        )
    )


@dataclass
class AbundanceNoExclusiveSpecies:
    description: str = "2 communities; neither contain exclusive species"
    counts: ndarray = field(
        default_factory=lambda: array(
            [
                [2, 4],
                [3, 1],
                [1, 5],
            ],
        )
    )
    subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 16, 4 / 16],
                [3 / 16, 1 / 16],
                [1 / 16, 5 / 16],
            ]
        )
    )
    metacommunity_abundance: ndarray = field(
        default_factory=lambda: array([[6 / 16], [4 / 16], [6 / 16]])
    )
    subcommunity_normalizing_constants: ndarray = field(
        default_factory=lambda: array([6 / 16, 10 / 16])
    )
    normalized_subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 6, 4 / 10],
                [3 / 6, 1 / 10],
                [1 / 6, 5 / 10],
            ]
        )
    )


@dataclass
class AbundanceMutuallyExclusive:
    description: str = "2 mutually exclusive communities"
    counts: ndarray = field(
        default_factory=lambda: array(
            [
                [2, 0],
                [3, 0],
                [0, 1],
                [0, 4],
            ],
        )
    )
    subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 10, 0 / 10],
                [3 / 10, 0 / 10],
                [0 / 10, 1 / 10],
                [0 / 10, 4 / 10],
            ]
        )
    )
    metacommunity_abundance: ndarray = field(
        default_factory=lambda: array([[2 / 10], [3 / 10], [1 / 10], [4 / 10]])
    )
    subcommunity_normalizing_constants: ndarray = field(
        default_factory=lambda: array([5 / 10, 5 / 10])
    )
    normalized_subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 5, 0 / 5],
                [3 / 5, 0 / 5],
                [0 / 5, 1 / 5],
                [0 / 5, 4 / 5],
            ]
        )
    )


@dataclass
class AbundanceOneSubcommunity:
    description: str = "one community"
    counts: ndarray = field(
        default_factory=lambda: array(
            [
                [2],
                [5],
                [3],
            ],
        )
    )
    subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 10],
                [5 / 10],
                [3 / 10],
            ]
        )
    )
    metacommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 10],
                [5 / 10],
                [3 / 10],
            ]
        )
    )
    subcommunity_normalizing_constants: ndarray = field(
        default_factory=lambda: array([10 / 10])
    )
    normalized_subcommunity_abundance: ndarray = field(
        default_factory=lambda: array(
            [
                [2 / 10],
                [5 / 10],
                [3 / 10],
            ]
        )
    )


@mark.parametrize(
    "counts, expected",
    [
        (counts_array_3by2(), AbundanceForDiversity),
        (counts_dataframe_3by2, AbundanceForDiversity),
    ],
)
def test_make_abundance(counts, expected):
    abundance = make_abundance(counts=counts)
    assert isinstance(abundance, expected)


def test_make_abundance_not_implemented():
    with raises(NotImplementedError):
        make_abundance(counts=1)


@mark.parametrize(
    "test_case",
    [
        AbundanceExclusiveSpecies(),
        AbundanceSingleExclusiveSpecies(),
        AbundanceNoExclusiveSpecies(),
        AbundanceMutuallyExclusive(),
        AbundanceOneSubcommunity(),
    ],
)
class TestAbundance:
    def test_make_subcommunity_abundance(self, test_case):
        abundance = make_abundance(counts=test_case.counts)
        assert allclose(
            abundance.subcommunity_abundance,
            test_case.subcommunity_abundance,
        )

    def test_metacommunity_abundance(self, test_case):
        abundance = make_abundance(counts=test_case.counts)
        assert allclose(
            abundance.metacommunity_abundance, test_case.metacommunity_abundance
        )

    def test_subcommunity_normalizing_constants(self, test_case):
        abundance = make_abundance(counts=test_case.counts)
        assert allclose(
            abundance.subcommunity_normalizing_constants,
            test_case.subcommunity_normalizing_constants,
        )

    def test_normalized_subcommunity_abundance(self, test_case):
        abundance = make_abundance(counts=test_case.counts)
        assert allclose(
            abundance.normalized_subcommunity_abundance,
            test_case.normalized_subcommunity_abundance,
        )
