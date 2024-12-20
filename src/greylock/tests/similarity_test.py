"""Tests for diversity.similarity"""

from collections import defaultdict
from numpy import (
    allclose,
    ndarray,
    array,
    dtype,
    memmap,
    inf,
    float32,
    zeros,
    identity,
    maximum,
    exp,
)
from numpy.linalg import norm
from pandas import DataFrame
import scipy.sparse
from pytest import fixture, raises, mark

from greylock.log import LOGGER
from greylock.similarity import (
    SimilarityIdentity,
    SimilarityFromArray,
    SimilarityFromDataFrame,
    SimilarityFromFile,
    IntersetSimilarityFromFile,
    SimilarityFromFunction,
    SimilarityFromSymmetricFunction,
    IntersetSimilarityFromFunction,
    weighted_similarity_chunk_nonsymmetric,
    weighted_similarity_chunk_symmetric,
)
from greylock import Metacommunity
from greylock.exceptions import InvalidArgumentError
from greylock.abundance import make_abundance


@fixture
def similarity_function():
    return lambda a, b: 1 / sum(a * b)


def similarity_from_distance(a, b):
    return exp(-1 * norm(a - b))


similarity_array_3by3_1 = array(
    [
        [1, 0.5, 0.1],
        [0.5, 1, 0.2],
        [0.1, 0.2, 1],
    ]
)
similarity_array_3by3_2 = array(
    [
        [1.0, 0.9, 0.1],
        [0.9, 1.0, 0.2],
        [0.1, 0.2, 1.0],
    ]
)
similarity_dataframe_3by3 = DataFrame(
    data=similarity_array_3by3_1,
    columns=["species_1", "species_2", "species_3"],
    index=["species_1", "species_2", "species_3"],
)
similarity_filecontent_3by3_tsv = (
    "species_1\tspecies_2\tspecies_3\n"
    "1.0\t0.5\t0.1\n"
    "0.5\t1.0\t0.2\n"
    "0.1\t0.2\t1.0\n"
)
similarities_filecontents_3by3_csv = (
    "species_1,species_2,species_3\n" "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
)
# fmt: off
similarities_filecontents_32by3_csv = (
    "species_1,species_2,species_3\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "1.0,0.5,0.1\n" "0.5,1.0,0.2\n" "0.1,0.2,1.0\n"
    "0.0,0.5,0.0\n" "0.5,0.0,0.2\n" "0.0,0.2,0.0\n"
    "0.0,0.5,0.0\n" "0.5,0.0,0.2\n" "0.0,0.2,0.0\n"
    "0.0,0.5,0.0\n" "0.5,0.0,0.2\n" "0.0,0.2,0.0\n"
    "0.0,0.5,0.0\n" "0.5,0.0,0.2\n" "0.0,0.2,0.0\n"
)
# fmt: on
similarity_sparse_entries = {
    "shape": (3, 3),
    "row": array(
        [
            0,
            1,
            2,
            0,
            1,
            0,
            2,
        ]
    ),
    "col": array(
        [
            0,
            1,
            2,
            1,
            0,
            2,
            0,
        ]
    ),
    "data": array(
        [
            1,
            1,
            1,
            0.5,
            0.5,
            0.3,
            0.3,
        ]
    ),
}
relative_abundance_3by1 = array([[1 / 1000], [1 / 10], [10]])
relative_abundance_3by2 = array([[1 / 1000, 1 / 100], [1 / 10, 1 / 1], [10, 100]])
relative_abundance_3by2_2 = array(
    [
        [0.7, 0.3],
        [0.1, 0.3],
        [0.2, 0.4],
    ]
)
weighted_abundances_3by1_1 = array([[1.051], [2.1005], [10.0201]])
weighted_abundances_3by1_2 = array([[0.35271989], [0.13459705], [0.0601738]])
weighted_abundances_3by1_3 = array([[0.35271989], [0.13459705], [0.0601738]])
weighted_abundances_3by2_1 = array(
    [[1.051, 10.51], [2.1005, 21.005], [10.0201, 100.201]]
)
weighted_abundances_3by2_2 = array(
    [
        [0.81, 0.61],
        [0.77, 0.65],
        [0.29, 0.49],
    ]
)
weighted_abundances_3by2_3 = (
    array(
        [
            [0.35271989, 3.52719894],
            [0.13459705, 1.34597047],
            [0.0601738, 0.60173802],
        ]
    ),
)
weighted_abundances_3by2_4 = array(
    [
        [1.46290476, 14.62904762],
        [0.48763492, 4.87634921],
        [0.20898639, 2.08986395],
    ]
)
X_3by1 = array([[1], [3], [7]])
X_3by2 = array([[1, 2], [3, 5], [7, 11]])


@fixture
def memmapped_similarity_matrix(tmp_path):
    memmapped = memmap(
        tmp_path / "similarity_matrix.npy",
        dtype=dtype("f8"),
        mode="w+",
        offset=0,
        shape=similarity_array_3by3_1.shape,
        order="C",
    )
    memmapped[:, :] = similarity_array_3by3_1
    return memmapped


@fixture
def make_similarity_from_file(tmp_path):
    def make(
        filename="similarity_matrix.tsv",
        filecontent=similarity_filecontent_3by3_tsv,
        chunk_size=1,
        similarity_class=SimilarityFromFile,
    ):
        filepath = tmp_path / filename
        with open(filepath, "w") as file:
            file.write(filecontent)
        return similarity_class(similarity_file_path=filepath, chunk_size=chunk_size)

    return make


@mark.parametrize(
    "relative_abundance, expected, kwargs",
    [
        (relative_abundance_3by2, weighted_abundances_3by2_1, {}),
        (relative_abundance_3by2, weighted_abundances_3by2_1, {"chunk_size": 2}),
        (
            relative_abundance_3by2,
            weighted_abundances_3by2_1,
            {
                "filename": "similarity_matrix.csv",
                "filecontent": similarities_filecontents_3by3_csv,
            },
        ),
        (relative_abundance_3by1, weighted_abundances_3by1_1, {}),
    ],
)
def test_weighted_abundances(
    relative_abundance, expected, kwargs, make_similarity_from_file
):
    similarity = make_similarity_from_file(**kwargs)
    assert allclose(similarity.weighted_abundances(relative_abundance), expected)


def test_nonsquare_from_file(make_similarity_from_file):
    sim = make_similarity_from_file(similarity_class=IntersetSimilarityFromFile)
    counts = array([[1], [1], [1]])
    with raises(InvalidArgumentError):
        Metacommunity(counts, sim).to_dataframe(viewpoint=0)


def test_interset_from_file(make_similarity_from_file):
    sim = make_similarity_from_file(
        filecontent=similarities_filecontents_32by3_csv,
        chunk_size=4,
        similarity_class=IntersetSimilarityFromFile,
    )
    counts = array([[50], [25], [25]])
    abundance_object = make_abundance(counts, False)
    # fmt: off
    expected = array([[0.65, 0.55, 0.35,
                      0.65, 0.55, 0.35,
                      0.65, 0.35,
                      0.65, 0.55, 0.35,
                      0.65, 0.55, 0.35,
                      0.65, 0.55, 0.35,
                      0.65, 0.55, 0.35,
                      0.125, 0.3, 0.05,
                      0.125, 0.3, 0.05,
                      0.125, 0.3, 0.05,
                       0.125, 0.3, 0.05]]).T
    # fmt: on
    result = sim @ abundance_object
    assert allclose(result, expected)


@mark.parametrize(
    "similarity, simclass, relative_abundance, expected",
    [
        (
            similarity_dataframe_3by3,
            SimilarityFromDataFrame,
            relative_abundance_3by2,
            weighted_abundances_3by2_1,
        ),
        (
            similarity_dataframe_3by3,
            SimilarityFromDataFrame,
            relative_abundance_3by1,
            weighted_abundances_3by1_1,
        ),
        (
            similarity_array_3by3_1,
            SimilarityFromArray,
            relative_abundance_3by1,
            weighted_abundances_3by1_1,
        ),
        (
            similarity_array_3by3_2,
            SimilarityFromArray,
            relative_abundance_3by2_2,
            weighted_abundances_3by2_2,
        ),
    ],
)
def test_weighted_abundances_from_array(
    similarity, simclass, relative_abundance, expected
):
    similarity = simclass(similarity=similarity)
    weighted_abundances = similarity.weighted_abundances(
        relative_abundance=relative_abundance
    )
    assert allclose(weighted_abundances, expected)


def test_weighted_abundances_from_memmap(memmapped_similarity_matrix):
    similarity = SimilarityFromArray(similarity=memmapped_similarity_matrix)
    weighted_abundances = similarity.weighted_abundances(
        relative_abundance=relative_abundance_3by2
    )
    assert allclose(weighted_abundances, weighted_abundances_3by2_1)


@mark.parametrize(
    "relative_abundance, X, chunk_size, expected",
    [
        (relative_abundance_3by2, X_3by2, 2, weighted_abundances_3by2_3),
        (relative_abundance_3by2, X_3by1, 1, weighted_abundances_3by2_4),
        (relative_abundance_3by1, X_3by2, 4, weighted_abundances_3by1_2),
        (relative_abundance_3by1, X_3by2, 2, weighted_abundances_3by1_2),
    ],
)
def test_weighted_abundances_from_function(
    relative_abundance, similarity_function, X, chunk_size, expected
):
    similarity = SimilarityFromFunction(
        func=similarity_function, X=X, chunk_size=chunk_size
    )
    weighted_abundances = similarity.weighted_abundances(
        relative_abundance=relative_abundance
    )
    assert allclose(weighted_abundances, expected)


def test_weighted_similarity_chunk(similarity_function):
    chunk_index, chunk = weighted_similarity_chunk_nonsymmetric(
        similarity=similarity_function,
        X=X_3by2,
        Y=None,
        relative_abundance=relative_abundance_3by2,
        chunk_size=3,
        chunk_index=0,
    )
    assert chunk_index == 0
    assert allclose(chunk, weighted_abundances_3by2_3)


def make_array(spec, array_class=zeros):
    a = array_class(spec["shape"], dtype=float32)
    for i, value in enumerate(spec["data"]):
        a[spec["row"][i], spec["col"][i]] = value
    return a


def compare_dense_sparse(counts, dense_similarity, sparse_similarity):
    viewpoints = [0, 1, 2, inf]
    measures = (
        "alpha",
        "rho",
        "beta",
        "gamma",
        "normalized_alpha",
        "normalized_rho",
        "normalized_beta",
    )
    meta_dense = Metacommunity(counts, similarity=dense_similarity)
    meta_dense_df = meta_dense.to_dataframe(viewpoint=viewpoints, measures=measures)
    meta_sparse = Metacommunity(
        counts, similarity=SimilarityFromArray(sparse_similarity)
    )
    meta_sparse_df = meta_sparse.to_dataframe(viewpoint=viewpoints, measures=measures)
    assert meta_dense_df.equals(meta_sparse_df)


@mark.parametrize(
    "sparse_class",
    [
        scipy.sparse.bsr_array,
        scipy.sparse.coo_array,
        scipy.sparse.csc_array,
        scipy.sparse.csr_array,
        scipy.sparse.bsr_matrix,
        scipy.sparse.coo_matrix,
        scipy.sparse.csc_matrix,
        scipy.sparse.csr_matrix,
    ],
)
def test_sparse_similarity(sparse_class):
    spec = similarity_sparse_entries
    dense_similarity = make_array(spec)
    sparse_similarity = sparse_class(
        (spec["data"], (spec["row"], spec["col"])), shape=spec["shape"]
    )
    counts = DataFrame({"Medford": [3, 2, 0], "Somerville": [1, 4, 0]})
    compare_dense_sparse(counts, dense_similarity, sparse_similarity)


@mark.parametrize("sparse_class", [scipy.sparse.dia_array, scipy.sparse.dia_matrix])
def test_diag_sparse(sparse_class):
    data = array([[0.5] * 4, [1] * 4, [0.5] * 4])
    offsets = array([-1, 0, 1])
    dense_similarity = array(
        [
            [1.0, 0.5, 0.0, 0.0],
            [0.5, 1.0, 0.5, 0.0],
            [0.0, 0.5, 1.0, 0.5],
            [0.0, 0.0, 0.5, 1.0],
        ]
    )
    sparse_similarity = sparse_class((data, offsets), shape=(4, 4))
    counts = DataFrame({"Cambridge": [5, 2, 0, 9], "Boston": [2, 3, 3, 2]})
    compare_dense_sparse(counts, dense_similarity, sparse_similarity)


@mark.parametrize(
    "sparse_class",
    [
        scipy.sparse.lil_array,
        scipy.sparse.lil_matrix,
        scipy.sparse.dok_array,
        scipy.sparse.dok_matrix,
    ],
)
def test_incremental_sparse(sparse_class):
    spec = similarity_sparse_entries
    dense_similarity = make_array(spec)
    sparse_similarity = make_array(spec, sparse_class)
    counts = DataFrame({"Arlington": [23, 12, 8], "Watertown": [15, 14, 19]})
    compare_dense_sparse(counts, dense_similarity, sparse_similarity)


def another_similarity_func(row_i, row_j):
    i = row_i[0] / row_i[1]
    j = row_j[0] / row_j[1]
    return min(1.0, max(0.0, 0.4 * (3 - abs(i - j))))


symmetric_example_X = array([[0, 88], [10, 10], [200, 100], [99, 33]])
symmetric_example_abundance = array([[1, 0], [0, 1], [1, 0], [0, 10]])


@mark.parametrize(
    "chunk_index, expected",
    [
        (0, [[0.4, 0.8], [1.6, 4.0], [0.4, 0.8], [0.0, 0.4]]),
        (2, [[0.0, 0.0], [0.0, 0.0], [0.0, 8.0], [0.8, 0.0]]),
    ],
)
def test_weighted_similarity_chunk_symmetric(chunk_index, expected):
    result = weighted_similarity_chunk_symmetric(
        another_similarity_func,
        symmetric_example_X,
        symmetric_example_abundance,
        2,
        chunk_index,
    )
    assert allclose(result, array(expected))


def test_symmetric_similarity():
    expected = array([[1.4, 0.8], [1.6, 5.0], [1.4, 8.8], [0.8, 10.4]])
    obj = SimilarityFromSymmetricFunction(
        func=another_similarity_func,
        X=symmetric_example_X,
        chunk_size=2,
    )
    result = obj.weighted_abundances(symmetric_example_abundance)
    assert allclose(result, expected)


animal_features = DataFrame(
    {
        "breathes": [
            "air",
            "air",
            "air",
            "air",
            "water",
            "air",
            "water",
            "air",
            "air",
            "air",
            "air",
            "air",
        ],
        "covering": [
            "fur",
            "fur",
            "fur",
            "scales",
            "scales",
            "fur",
            "scales",
            "feathers",
            "fur",
            "fur",
            "scales",
            "fur",
        ],
        "diet": [
            "meat",
            "omni",
            "plants",
            "meat",
            "omni",
            "plants",
            "omni",
            "omni",
            "meat",
            "plants",
            "plants",
            "plants",
        ],
        "n_legs": [
            4,
            4,
            4,
            0,
            0,
            4,
            0,
            2,
            4,
            4,
            4,
            2,
        ],
    },
    index=[
        "cat",
        "dog",
        "rabbit",
        "snake",
        "goldfish",
        "gerbil",
        "betafish",
        "chicken",
        "tiger",
        "giraffe",
        "turtle",
        "monkey",
    ],
)

animal_communities = DataFrame(
    {
        "zoo": [
            0,
            0,
            2,
            6,
            0,
            0,
            0,
            0,
            2,
            2,
            4,
            5,
        ],
        "petco": [
            4,
            6,
            4,
            4,
            40,
            8,
            30,
            0,
            0,
            0,
            1,
            0,
        ],
        "shelter": [
            10,
            15,
            1,
            1,
            1,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
        ],
    },
    index=[
        "cat",
        "dog",
        "rabbit",
        "snake",
        "goldfish",
        "gerbil",
        "betafish",
        "chicken",
        "tiger",
        "giraffe",
        "turtle",
        "monkey",
    ],
)


def feature_similarity(animal_i, animal_j):
    if animal_i.breathes != animal_j.breathes:
        return 0.0
    if animal_i.covering == animal_j.covering:
        if animal_i.Index == animal_j.Index:
            result = 1
        else:
            result = 0.98
    else:
        result = 0.5
    if animal_i.n_legs != animal_j.n_legs:
        result *= 0.5
    if animal_i.diet != animal_j.diet:
        if "omni" in (animal_i.diet, animal_j.diet):
            result *= 0.9
        else:
            result *= 0.7
    return result


def animal_similarity_matrix():
    columns = defaultdict(list)
    index = []
    results = []
    for i, animal_i in enumerate(animal_features.itertuples()):
        index.append(animal_i.Index)
        for j, animal_j in enumerate(animal_features.itertuples()):
            s = feature_similarity(animal_i, animal_j)
            columns[animal_j.Index].append(s)
            results.append((s, animal_i.Index, animal_j.Index))
    return DataFrame(columns, index=index)


def test_feature_similarity():
    measures = [
        "alpha",
        "rho",
        "beta",
        "gamma",
        "normalized_alpha",
        "normalized_rho",
        "normalized_beta",
        "rho_hat",
    ]
    viewpoints = [0, 1, 2, inf]
    m = Metacommunity(
        animal_communities,
        similarity=SimilarityFromDataFrame(animal_similarity_matrix()),
    )
    df1 = m.to_dataframe(viewpoint=viewpoints, measures=measures).set_index(
        ["community", "viewpoint"]
    )
    m = Metacommunity(
        animal_communities,
        SimilarityFromFunction(
            func=feature_similarity,
            X=animal_features,
            chunk_size=4,
        ),
    )
    df2 = m.to_dataframe(viewpoint=viewpoints, measures=measures).set_index(
        ["community", "viewpoint"]
    )
    assert allclose(df1.to_numpy(), df2.to_numpy())
    m = Metacommunity(
        animal_communities,
        SimilarityFromSymmetricFunction(
            func=feature_similarity,
            X=animal_features,
            chunk_size=4,
        ),
    )
    df3 = m.to_dataframe(viewpoint=viewpoints, measures=measures).set_index(
        ["community", "viewpoint"]
    )
    assert allclose(df1.to_numpy(), df3.to_numpy())


def test_nonsymmetric():
    """
    SimilarityFromFunction will yield, if the function is not
    actually symmetric, a matrix that is not a reflection of itself
    across the diagonal.
    """

    def nonsym_similarity_function(species_i, species_j):
        diff = maximum((species_i - species_j), 0)
        return 1 / (1 + norm(diff) / 100)

    X = array(
        [[54, 200, 45, 123], [55, 67, 44, 99], [25, 145, 56, 12], [154, 98, 55, 98]]
    )
    counts = identity(4)
    sim = SimilarityFromFunction(nonsym_similarity_function, X)
    matrix = zeros(shape=(4, 4))
    for i in range(4):
        for j in range(4):
            matrix[i, j] = nonsym_similarity_function(X[i], X[j])
    result = sim.weighted_abundances(counts)
    assert allclose(result, matrix)
    for i in range(4):
        for j in range(i + 1, 4):
            assert result[i, j] != result[j, i]


@fixture
def callcounter():
    return defaultdict(int)


@mark.parametrize(
    "sim, key, expected_count",
    [
        [SimilarityIdentity(), "identity", 3],
        ["array", "array", 3],
        ["df", "df", 3],
        ["file", "file", 1],
        [
            SimilarityFromFunction(
                similarity_from_distance, X=array([[1, 2], [2, 1], [3, 0]])
            ),
            "func",
            1,
        ],
        [
            SimilarityFromSymmetricFunction(
                similarity_from_distance, X=array([[1, 2], [2, 1], [3, 0]])
            ),
            "sfunc",
            1,
        ],
    ],
)
def test_compuation_count(
    sim, key, expected_count, callcounter, make_similarity_from_file
):
    """
    Test that we unify the abundance array and call weighted_abundances only once
    for the subclasses of Similarity that are expensive to run, but don't bother
    with that for trivial subclasses.
    """

    def count_decorator(f, counter, key):
        def wrapper(*args, **kwds):
            counter[key] += 1
            return f(*args, **kwds)

        return wrapper

    abundances = array([[1, 2], [2, 5], [9, 3]])
    if sim == "file":
        sim = make_similarity_from_file()
    elif sim == "array":
        sim = SimilarityFromArray(similarity_array_3by3_1)
    elif sim == "df":
        sim = SimilarityFromDataFrame(similarity_dataframe_3by3)

    sim.weighted_abundances = count_decorator(sim.weighted_abundances, callcounter, key)
    m = Metacommunity(abundances, sim)
    m.metacommunity_diversity(viewpoint=1, measure="alpha")
    m.to_dataframe(viewpoint=0)

    assert callcounter[key] == expected_count


@mark.parametrize(
    "X, Y, abundance, expected",
    [
        [
            array([[2, 1, 5], [4, 3, 2]]),
            array([[2, 0, 5], [4, 2, 1], [2, 0, 5]]),
            array([[0.5, 0.1], [0.25, 0.1], [0.25, 0.8]]),
            array([[0.27846671, 0.33211435], [0.06766633, 0.03257625]]),
        ],
        [
            array([[1, 0, 0]]),
            array([[1, 0, 0], [0, 1, 0], [0, 1, 1]]),
            array([[1.0, 0.8, 0.0], [0.0, 0.1, 0.5], [0.0, 0.1, 0.5]]),
            array([[1.0, 0.84200379, 0.21001897]]),
        ],
    ],
)
def test_interset_similarity(X, Y, abundance, expected):
    sim = IntersetSimilarityFromFunction(similarity_from_distance, X, Y)
    result = sim.weighted_abundances(abundance)
    assert allclose(result, expected)
    abundance_obj = make_abundance(abundance, False)
    result = sim @ abundance_obj
    assert allclose(result, expected)


def test_square_similarity_requirement():
    sim = SimilarityFromArray(array([[1.0, 0.4, 0.6], [0.2, 1.0, 0.9]]))
    counts = array([[1], [1], [1]])
    with raises(InvalidArgumentError):
        Metacommunity(counts, sim).to_dataframe(viewpoint=0)


def test_interset_diversity_forbidden():
    sim = IntersetSimilarityFromFunction(
        similarity_from_distance,
        X=array([[1, 0, 0]]),
        Y=array([[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]]),
    )
    counts = array([[1, 1, 1, 1, 1]])
    with raises(InvalidArgumentError):
        Metacommunity(counts, sim).to_dataframe(viewpoint=0)


def test_matmul():
    similarity_df = DataFrame(
        {"bee": [0.6, 0.22], "butterfly": [0.55, 0.27], "lobster": [0.45, 0.28]},
        index=["ladybug", "fish"],
    )
    similarity = SimilarityFromDataFrame(similarity_df)
    abundance = make_abundance(array([[5000], [2000], [3000]]), False)
    expected = array([[0.545], [0.248]])
    relative_abundance = similarity @ abundance
    assert allclose(expected, relative_abundance)
