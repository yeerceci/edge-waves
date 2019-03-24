% function to put all the elements of an MxN matrix, X, onto the diagonal of an MNxMN matrix, Y.

function Y = longdiag(X)

Y = diag(X(:));