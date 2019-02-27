# coding=utf-8
# Copyright 2019 The Tensor2Tensor Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
python s1_build_token_vocab.py datadir
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import os
import sys
from tensor2tensor.data_generators import text_encoder
import tensorflow as tf


EOS = text_encoder.EOS
data_dir=sys.argv[1]
# data_dir="data_token_sub"


def _read_words(filename):
  """Reads words from a file."""
  with tf.gfile.GFile(filename, "r") as f:
    if sys.version_info[0] >= 3:
      return f.read().replace("\n", " %s " % EOS).split()
    else:
      return f.read().decode("utf-8").replace("\n", " %s " % EOS).split()


def _build_vocab(filenames, vocab_path):
  """Reads a file to build a vocabulary of `vocab_size` most common words.
   The vocabulary is sorted by occurrence count and has one word per line.
   Originally from:
   https://github.com/tensorflow/models/blob/master/tutorials/rnn/ptb/reader.py
  Args:
    filename: file to read list of words from.
    vocab_path: path where to save the vocabulary.
    vocab_size: size of the vocabulary to generate.
  """
  data=_read_words(filenames[0])
  for file in filenames[1:]:
    data=data+_read_words(file)

  counter = collections.Counter(data)
  count_pairs = sorted(counter.items(), key=lambda x: (-x[1], x[0]))
  words, _ = list(zip(*count_pairs))
  # words = words[:vocab_size]
  a=list(words)
  a.remove('<EOS>')
  with open(vocab_path, "w") as f:
    f.write("<UNK>\n")
    f.write("<EOS>\n")
    f.write("\n".join(a))

files=[]
files.append(os.path.join(data_dir, "train_sources"))
files.append(os.path.join(data_dir, "train_targets"))
if os.path.exists(os.path.join(data_dir, "test_sources")):
  files.append(os.path.join(data_dir, "test_sources"))
if os.path.exists(os.path.join(data_dir, "test_r1_sources")):
  files.append(os.path.join(data_dir, "test_r1_sources"))

vocab_file = os.path.join(data_dir, "vocab.token")
_build_vocab(files, vocab_file)









