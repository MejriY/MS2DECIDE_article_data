import pytest
import pandas as pd

def test_read_tsv():
    # Test reading the file in resource
    path = "tests/resources/75b94088c0f24a0eb8b064b380300649.tsv"
    df = pd.read_csv(path, sep="\t")
    assert df.shape == (1, 2)
