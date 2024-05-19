import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

## filter noise from data

def pca_on_sim_data(dataset_name="data4pca.mat"):

    data = loadmat(f"data/{dataset_name}")
    # X must have shape (num_dimensions, num_observations)
    X = data['highdim_dat'].T
    num_dims = X.shape[0]
    num_samples = X.shape[1]

    Z = X - np.mean(X, axis=1).reshape((num_dims, -1))
    sigma = np.cov(Z, rowvar=True)

    # get eigenvectors & eigenvalues of the covariance matrix
    eigenvalues = np.linalg.eig(sigma)[0]
    eigenvectors = np.linalg.eig(sigma)[1]

    # sort from high to low
    index = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[index]
    print(f"eigenvalues: {eigenvalues}")
    principal_components = eigenvectors[:, index]
    print(f"PCs: {principal_components}")

    lambda_sum = np.sum(eigenvalues)
    PC1_PC2_variance = (eigenvalues[0] + eigenvalues[1]) / lambda_sum
    print(f"PC1 and PC2 total variance explained: {PC1_PC2_variance}")
    Z_in_PC_basis = principal_components.T @ Z

    Z_in_PC_basis_new = principal_components[:, 0:2].T @ Z

    plt.scatter(Z_in_PC_basis_new[0, :], Z_in_PC_basis_new[1, :])
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Z in the basis of {P1, PC2}")
    plt.show()

pca_on_sim_data()