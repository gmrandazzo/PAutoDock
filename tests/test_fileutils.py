#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tests_fileutils.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

unit tests for fileutils.py

"""
import pytest
from unittest.mock import patch, mock_open
import requests
import tarfile
from pautodock.fileutils import (
    get_bin_path,
    download_file,
    extract_tar_gz
)


def test_get_bin_path_linux_first_path():
    with patch('platform.system', return_value='Linux'):
        with patch('pathlib.Path.exists', side_effect=[True, False]):
            assert get_bin_path('testbin') == '/usr/bin/'

def test_get_bin_path_linux_second_path():
    with patch('platform.system', return_value='Linux'):
        with patch('pathlib.Path.exists', side_effect=[False, True]):
            assert get_bin_path('testbin') == '/usr/local/bin'

def test_get_bin_path_darwin():
    with patch('platform.system', return_value='Darwin'):
        with patch('pathlib.Path.exists', return_value=True):
            assert get_bin_path('testbin') == '/opt/homebrew/bin/'

def test_get_bin_path_unsupported_platform():
    with patch('platform.system', return_value='Windows'):
        with pytest.raises(ValueError, match="Platform not supported"):
            get_bin_path('testbin')

def test_get_bin_path_not_found_linux():
    with patch('platform.system', return_value='Linux'):
        with patch('pathlib.Path.exists', return_value=False):
            with pytest.raises(ValueError, match="Unable to find testbin installed."):
                get_bin_path('testbin')

def test_get_bin_path_not_found_darwin():
    with patch('platform.system', return_value='Darwin'):
        with patch('pathlib.Path.exists', return_value=False):
            with pytest.raises(ValueError, match="Unable to find testbin installed."):
                get_bin_path('testbin')

@pytest.mark.parametrize("bin_name", ["obabel", "vina", "custom_bin"])
def test_get_bin_path_different_binaries(bin_name):
    with patch('platform.system', return_value='Linux'):
        with patch('pathlib.Path.exists', return_value=True):
            assert get_bin_path(bin_name) == '/usr/bin/'

def test_get_bin_path_case_sensitivity():
    with patch('platform.system', return_value='Linux'):
        with patch('pathlib.Path.exists', side_effect=[False, True]):
            assert get_bin_path('TestBin') == '/usr/local/bin'

# Test download_file function
@patch('requests.get')
@patch('builtins.open', new_callable=mock_open)
def test_download_file_success(mock_file, mock_get):
    # Setup
    mock_response = mock_get.return_value
    mock_response.iter_content.return_value = [b'chunk1', b'chunk2']
    mock_response.raise_for_status.return_value = None

    # Execute
    result = download_file('http://example.com/file', 'local_file')

    # Assert
    assert result is True
    mock_get.assert_called_once_with('http://example.com/file', stream=True)
    mock_file.assert_called_once_with('local_file', 'wb')
    mock_file().write.assert_any_call(b'chunk1')
    mock_file().write.assert_any_call(b'chunk2')

@patch('requests.get')
def test_download_file_failure(mock_get):
    # Setup
    mock_get.side_effect = requests.RequestException('Network error')

    # Execute
    result = download_file('http://example.com/file', 'local_file')

    # Assert
    assert result is False
    mock_get.assert_called_once_with('http://example.com/file', stream=True)

# Test extract_tar_gz function
@patch('tarfile.open')
def test_extract_tar_gz_success(mock_tarfile):
    # Setup
    mock_tar = mock_tarfile.return_value.__enter__.return_value

    # Execute
    result = extract_tar_gz('file.tar.gz', 'extract_dir')

    # Assert
    assert result is True
    mock_tarfile.assert_called_once_with('file.tar.gz', 'r:gz')
    mock_tar.extractall.assert_called_once_with(path='extract_dir')

@patch('tarfile.open')
def test_extract_tar_gz_failure(mock_tarfile):
    # Setup
    mock_tarfile.side_effect = tarfile.TarError('Extraction error')

    # Execute
    result = extract_tar_gz('file.tar.gz', 'extract_dir')

    # Assert
    assert result is False
    mock_tarfile.assert_called_once_with('file.tar.gz', 'r:gz')