import numpy as np

fullcov = np.zeros(shape=(90,90))

partcov = np.array([ \
(1., 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99), \
(0.99, 1., 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99), \
(0.99, 0.99, 1., 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99), \
(0.99, 0.99, 0.99, 1., 0.99, 0.99, 0.99, 0.99, 0.99, 0.99), \
(0.99, 0.99, 0.99, 0.99, 1., 0.99, 0.99, 0.99, 0.99, 0.99), \
(0.99, 0.99, 0.99, 0.99, 0.99, 1., 0.99, 0.99, 0.99, 0.99), \
(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1., 0.99, 0.99, 0.99), \
(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1., 0.99, 0.99), \
(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1., 0.99), \
(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1.) ])

L_matrix = np.identity(9)

for i in range(9):
    for j in range(9):
        fullcov[i*10:(i+1)*10,j*10:(j+1)*10] = partcov * L_matrix[i,j]

C = np.linalg.cholesky(fullcov)
a = np.random.randn(90)

#print(np.round(fullcov,3))
print(np.round(C,3))
#print(np.dot(a,C))
