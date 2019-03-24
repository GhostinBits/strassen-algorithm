#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

void submatrix(double * result, const double * original, int y, int x,  // Creator of submatrix according to its original's offsets
               int offset_to_ystart, int offset_to_xstart, int offset_to_yend, int offset_to_xend) {
    int suby = y - offset_to_ystart - offset_to_yend;
    int subx = x - offset_to_xstart - offset_to_xend;
    for (int i = 0, i_w_o = offset_to_ystart; i < suby; i++) {
        for (int j = 0, j_w_o = offset_to_xstart; j < subx; j++) {
            result[i * subx + j] = original[i_w_o * x + j_w_o];
        }
    }
}

void add(double * result, const double * matrix_1, int y, int x, const double * matrix_2) {
    for (int i = 0; i < y; i++) {
        for (int j = 0; j < x; j++) {
            result[i * x + j] = matrix_1[i * x + j] + matrix_2[i * x + j];
        }
    }
}

void sub(double * result, const double * matrix_1, int y, int x, const double * matrix_2) {
    for (int i = 0; i < y; i++) {
        for (int j = 0; j < x; j++) {
            result[i * x + j] = matrix_1[i * x + j] - matrix_2[i * x + j];
        }
    }
}


void mul(double * result, const double * matrix_1, int y_1, int x_1, const double * matrix_2, int y_2, int x_2) {
    for (int i = 0; i < y_1; i++) {
        for (int j = 0; j < x_2; j++) {
            result[i * x_2 + j] = 0;
            for (int k = 0; k < x_1; k++) {
                result[i * x_2 + j] += matrix_1[i * x_1 + k] * matrix_2[k * x_2 + j];
            }
        }
    }
}

void print_matrix(double * matrix, int y, int x) {
    for (int i = 0; i < y; i++) {
        for (int j = 0; j < x; j++) {
            cout << matrix[i * x + j] << " ";
        }
        cout << endl;
    }
}


void strassen(double * result, const double * matrix_1, int y_1, int x_1, const double * matrix_2, int y_2, int x_2) {
    double A11[(y_1 / 2 * (x_1 / 2))];
    double A12[(y_1 / 2 * (x_1 / 2))];
    double A21[(y_1 / 2 * (x_1 / 2))];
    double A22[(y_1 / 2 * (x_1 / 2))];
    double B11[(y_2 / 2 * (x_2 / 2))];
    double B12[(y_2 / 2 * (x_2 / 2))];
    double B21[(y_2 / 2 * (x_2 / 2))];
    double B22[(y_2 / 2 * (x_2 / 2))];

    double S1[(y_2 / 2) * (x_2 / 2)];
    double S2[(y_1 / 2) * (x_1 / 2)];
    double S3[(y_1 / 2) * (x_1 / 2)];
    double S4[(y_2 / 2) * (x_2 / 2)];
    double S5[(y_1 / 2) * (x_1 / 2)];
    double S6[(y_2 / 2) * (x_2 / 2)];
    double S7[(y_1 / 2) * (x_1 / 2)];
    double S8[(y_2 / 2) * (x_2 / 2)];
    double S9[(y_1 / 2) * (x_1 / 2)];
    double S10[(y_2 / 2) * (x_2 / 2)];

    double P1[(y_1 / 2) * (x_2 / 2)];
    double P2[(y_1 / 2) * (x_2 / 2)];
    double P3[(y_1 / 2) * (x_2 / 2)];
    double P4[(y_1 / 2) * (x_2 / 2)];
    double P5[(y_1 / 2) * (x_2 / 2)];
    double P6[(y_1 / 2) * (x_2 / 2)];
    double P7[(y_1 / 2) * (x_2 / 2)];

    double Q1[(y_1 / 2) * (x_2 / 2)];
    double Q2[(y_1 / 2) * (x_2 / 2)];

    double C11[(y_1 / 2) * (x_2 / 2)];
    double C12[(y_1 / 2) * (x_2 / 2)];
    double C21[(y_1 / 2) * (x_2 / 2)];
    double C22[(y_1 / 2) * (x_2 / 2)];

    submatrix(A11, matrix_1, y_1, x_1, 0, 0, y_1 / 2, x_1 / 2);
    submatrix(A12, matrix_1, y_1, x_1, 0, x_1 / 2, y_1 / 2, 0);
    submatrix(A21, matrix_1, y_1, x_1, y_1 / 2, 0, 0, x_1 / 2);
    submatrix(A22, matrix_1, y_1, x_1, y_1 / 2, x_1 / 2, 0, 0);
    submatrix(B11, matrix_2, y_2, x_2, 0, 0, y_2 / 2, x_2 / 2);
    submatrix(B12, matrix_2, y_2, x_2, 0, x_2 / 2, y_2 / 2, 0);
    submatrix(B21, matrix_2, y_2, x_2, y_2 / 2, 0, 0, x_2 / 2);
    submatrix(B22, matrix_2, y_2, x_2, y_2 / 2, x_2 / 2, 0, 0);

    sub(S1, B12, y_2 / 2, x_2 / 2, B22);
    add(S2, A11, y_1 / 2, x_1 / 2, A12);
    add(S3, A21, y_1 / 2, x_1 / 2, A22);
    sub(S4, B21, y_2 / 2, x_2 / 2, B11);
    add(S5, A11, y_1 / 2, x_1 / 2, A22);
    add(S6, B11, y_2 / 2, x_2 / 2, B22);
    sub(S7, A12, y_1 / 2, x_1 / 2, A22);
    add(S8, B21, y_2 / 2, x_2 / 2, B22);
    sub(S9, A11, y_1 / 2, x_1 / 2, A21);
    add(S10, B11, y_2 / 2, x_2 / 2, B12);

    mul(P1, A11, y_1 / 2, x_1 / 2, S1, y_2 / 2, x_2 / 2);
    mul(P2, S2, y_1 / 2, x_1 / 2, B22, y_2 / 2, x_2 / 2);
    mul(P3, S3, y_1 / 2, x_1 / 2, B11, y_2 / 2, x_2 / 2);
    mul(P4, A22, y_1 / 2, x_1 / 2, S4, y_2 / 2, x_2 / 2);
    mul(P5, S5, y_1 / 2, x_1 / 2, S6, y_2 / 2, x_2 / 2);
    mul(P6, S7, y_1 / 2, x_1 / 2, S8, y_2 / 2, x_2 / 2);
    mul(P7, S9, y_1 / 2, x_1 / 2, S10, y_2 / 2, x_2 / 2);

    add(Q1, P5, y_1 / 2, x_2 / 2, P4);
    sub(Q2, P2, y_1 / 2, x_2 / 2, P6);
    sub(C11, Q1, y_1 / 2, x_2 / 2, Q2);

    add(C12, P1, y_1 / 2, x_2 / 2, P2);

    add(C21, P3, y_1 / 2, x_2 / 2, P4);

    add(Q1, P5, y_1 / 2, x_2 / 2, P1);
    add(Q2, P3, y_1 / 2, x_2 / 2, P7);
    sub(C22, Q1, y_1 / 2, x_2 / 2, Q2);

    for (int i = 0; i < y_1 / 2; i++) {
        for (int j = 0; j < x_2 / 2; j++) {
            result[i * x_2 + j] = C11[i * x_2 + j];
        }
    }
    for (int i = 0; i < y_1 / 2; i++) {
        for (int j = x_2 / 2, j_wo_offset = 0; j < x_2; j++, j_wo_offset++) {
            result[i * x_2 + j] = C12[i * x_2 + j_wo_offset];
        }
    }
    for (int i = y_1 / 2, i_wo_offset = 0; i < y_1; i++, i_wo_offset++) {
        for (int j = 0; j < x_2 / 2; j++) {
            result[i * x_2 + j] = C21[i_wo_offset * x_2 + j];
        }
    }
    for (int i = y_1 / 2, i_wo_offset = 0; i < y_1; i++, i_wo_offset++) {
        for (int j = x_2 / 2, j_wo_offset = 0; j < x_2; j++, j_wo_offset++) {
            result[i * x_2 + j] = C22[i_wo_offset * x_2 + j_wo_offset];
        }
    }
}



int main() {
    double A[100];
    double B[100];
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            A[i * 10 + j] = 1;
            B[i * 10 + j] = 1;
        }
    }
    auto * C = new double[100];
    strassen(C, A, 10, 10, B, 10, 10);
    print_matrix(C, 10, 10);

    double D[400];
    double E[400];
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            D[i * 10 + j] = M_PI;
            E[i * 10 + j] = exp(1);
        }
    }
    auto * F = new double[400];
    strassen(F, D, 20, 20, E, 20, 20);
    print_matrix(F, 20, 20);

    return 0;
}
