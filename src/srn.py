# A simple recurrent neural network.

# import stuff
import os
import numpy as np
import sys
from keras.layers import Dense, Activation
from keras.layers.recurrent import SimpleRNN
from keras.models import Sequential

hidden_neurons = 42
optimizer = 'sgd'
batch_size = 50
error_function = 'mean_squared_error'
output_nonlinearity = 'softmax'
cycles = 5
epochs_per_cycle = 3


def find_files(path, tag):
  file_paths = []
  for root, dirs, files in os.walk(path):
    for file in files:
      if tag in file:
        path = os.path.join(root,file)
        file_paths.append(path)
  return file_paths


file_path = find_files('../datasets', 'VDJ')[0]

cdr3s = ''
contents = open(file_path).read().splitlines()
for content in contents[1:]:
  parts = content.split('\t')
  cdr3 = parts[1]
  cdr3s += cdr3

cdr3_list = list(cdr3s)
chars = sorted(set(cdr3_list))
n_chars = len(chars)
char2index = dict((char,i) for i,char in enumerate(chars))
index2char = dict((i,char) for i,char in enumerate(chars))


def target_input_pairs(corpus, context):
  '''
  creates target input pairs from corpus using the given context
  '''
  inputs = []
  targets = []
  for i in range(len(corpus)-context):
    input = corpus[i:i+context]
    target = corpus[i+context]
    inputs.append(input)
    targets.append(target)
  return inputs, targets

context = 6
inputs,targets = target_input_pairs(cdr3_list,context)

input_vectors = np.zeros((len(inputs), context, n_chars ),dtype=np.int16)
target_vectors=  np.zeros((len(inputs), n_chars), dtype=np.int16)

for i, input in enumerate(inputs[:2]):
  for j, char in enumerate(input):
    input_vectors[i,j,char2index[char]] = 1
    tchar = targets[i]
    target_vectors[i,char2index[tchar]] = 1

model = Sequential()
model.add(SimpleRNN(hidden_neurons, input_shape=(context,n_chars), unroll = True))
model.add(Dense(n_chars))
model.add(Activation(output_nonlinearity))
model.compile(loss=error_function, optimizer=optimizer)
print(model.weights)
for cycle in range(cycles):
  print('>-<' * 50)
  print('cycle %d' % cycle)
  #model.fit(input_vectors, target_vectors, batch_size=batch_size, epochs=epochs_per_cycle)
  test_index = np.random.randint(len(inputs))
  test_chars = inputs[test_index]
  print('generating from sample %s, %s' % (test_index, test_chars))
  input_for_test = np.zeros((1,context,n_chars))
  for i,w in enumerate(test_chars):
    input_for_test[0,i,char2index[w]] = 1
  predictions_all_matrix = model.predict(input_for_test)[0]
  predicted_index = np.argmax(predictions_all_matrix)
  predicted_char = index2char[predicted_index]
  print(predicted_char)
 
