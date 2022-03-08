import numpy as np
import numpy
import matplotlib.pyplot as plt
import torch
from sklearn.datasets import make_blobs

x = torch.ones(1, requires_grad=True)
y = x + 2
z = y * y * 2
z.backward()

class Perceptron(torch.nn.Module):
    def __init__(self):
        super(Perceptron, self).__init__()
        self.fc = nn.Linear(1,1)
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        output = self.fc(x)
        output = self.relu(x) # instead of Heaviside step fn
        return output



class Feedforward(torch.nn.Module):
     def __init__(self, input_size, hidden_size):
         super(Feedforward, self).__init__()
         self.input_size = input_size
         self.hidden_size  = hidden_size
         self.fc1 = torch.nn.Linear(self.input_size, self.hidden_size)
         self.relu = torch.nn.ReLU()
         self.fc2 = torch.nn.Linear(self.hidden_size, 1)
         self.sigmoid = torch.nn.Sigmoid()
     def forward(self, x):
         hidden = self.fc1(x)
         relu = self.relu(hidden)
         output = self.fc2(relu)
         output = self.sigmoid(output)
         return output

# CREATE RANDOM DATA POINTS
def blob_label(y, label, loc): # assign labels
    target = numpy.copy(y)
    for l in loc:
        target[y == l] = label
    return target


x_train, y_train = make_blobs(n_samples=40, n_features=2, cluster_std=1.5, shuffle=True)
x_train = torch.FloatTensor(x_train)
y_train = torch.FloatTensor(blob_label(y_train, 0, [0]))
y_train = torch.FloatTensor(blob_label(y_train, 1, [1,2,3]))
x_test, y_test = make_blobs(n_samples=10, n_features=2, cluster_std=1.5, shuffle=True)
x_test = torch.FloatTensor(x_test)
y_test = torch.FloatTensor(blob_label(y_test, 0, [0]))
y_test = torch.FloatTensor(blob_label(y_test, 1, [1,2,3]))

model = Feedforward(2, 10)
criterion = torch.nn.BCELoss()
optimizer = torch.optim.SGD(model.parameters(), lr = 0.01)

model.eval()
y_pred = model(x_test)
before_train = criterion(y_pred.squeeze(), y_test)
print('Test loss before training' , before_train.item())






