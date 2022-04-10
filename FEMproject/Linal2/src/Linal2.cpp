
#include "Linal2.h"

using namespace std;


MyArray::MyArray() {
	this->array_size = 0;
	this->p = nullptr;
}

MyArray::MyArray(int array_size) {
    if (array_size == 0)
        MyArray();
    else {
        this->array_size = array_size;
        this->p = new float[array_size];
        for (int i = 0; i < array_size; i++) {
            p[i] = 0.0;
        }
    }
}


MyArray::MyArray(const MyArray &a) {
	this->array_size = a.array_size;
	this->p = new float[array_size];
	for (int i = 0; i < array_size; i++) {
		p[i] = a.p[i];
	}
}

MyArray::MyArray(int array_size, float a, float b) {
    if (array_size == 0)
        MyArray();
    else {
        this->array_size = array_size;
        this->p = new float[array_size];

        float h = (b - a) / (array_size - 1);
        p[0] = a;
        for (int i = 1; i < array_size; ++i) {
            p[i] = p[i - 1] + h;
        }
    }
}


MyArray::~MyArray() {
	delete [] p;
}

void MyArray::Resize(int new_size) {
    if (array_size == 0) {
        this->array_size = new_size;
        this->p = new float[array_size];
        for (int i = 0; i < array_size; i++) {
            p[i] = 0.0;
        }
    }
}

void MyArray::Show() {
	for (int i = 0; i < array_size; i++) {
		cout << p[i] << " ";
	}
	cout << endl;
}

void MyArray::Set(int index, float value) {
	if (p != nullptr) {
		if((index >= 0) && (index < array_size)) {
			p[index] = value;
		}
	}
}

void MyArray::add(MyArray a) {
    assert(this->array_size != a.get_size());
    for (int i = 0; i < array_size; ++i) {
        p[i] += a[i];
    }
}

void MyArray::add_weighted(MyArray a, float v1, float v2) {
    assert(this->array_size != a.get_size());
    for (int i = 0; i < this->array_size; ++i) {
        p[i] = v1 * p[i] + v2 * a[i];
    }
}

void MyArray::scale(float value) {
    for (int i = 0; i < array_size; i++) {
        p[i] *= value;
    }
}

int MyArray::get_size() {
	return array_size;
}

float* MyArray::get_data() {
    return p;
}

void MyArray::zap() {
	for (int i = 0; i < array_size; i++) {
		p[i] = rand() % 10;
	}
}
	
float & MyArray::operator [](int index) {
	return 	p[index];
}
	
MyArray MyArray::operator =(const MyArray &a) {
	this->array_size = a.array_size;
	this->p = new float[array_size];

   	for (int i = 0; i < array_size; i++) {
      		p[i] = a.p[i];
	}

	return *this;
}


float MyArray::norma() {
	float res = 0;
	for (int i = 0; i < array_size; i++) {
		res += p[i] * p[i];	
		//cout<<"res= "<<res<<endl;
	}
	return sqrt(abs(res));
}

void MyArray::grad(float (*f)(MyArray x), MyArray &res, float eps) {
    float tau = 0.1 * sqrt(eps);
    //float fx = f(x);
    //cout<<"fx= "<<fx<<endl;
    for (int i = 0; i < array_size; i++) {
        float tmp1 = p[i];
        p[i] += tau;
		float fx1 = f(*this);
		p[i] = tmp1;	
	
		float tmp2 = p[i];
		p[i]-=tau;
		float fx2 = f(*this);
		p[i] = tmp2;

        res[i] = (fx1 - fx2) / (2 * tau);
    }
}


float MyArray::dot_product(MyArray v2) {
    assert(this->array_size == v2.get_size());

    float res = 0.0;
    for (int i = 0; i < this->array_size; ++i) {
        res += this->p[i] * v2[i];
    }

    return res;
}

MyArray MyArray::multiplyElementwise(MyArray a) {
    assert(this->array_size == a.get_size());

    MyArray res(this->array_size);

    for (int i = 0; i < this->array_size; ++i) {
        res[i] = this->p[i] * a[i];
    }

    return res;
}

MyArray MyArray::divideByElementwise(MyArray a) {
    assert(this->array_size == a.get_size());

    MyArray res(this->array_size);

    for (int i = 0; i < this->array_size; ++i) {
        res[i] = this->p[i] / a[i];
    }

    return res;
}

//void MyArray::multiplyElementwise(MyArray a) {
//    assert(this->array_size == a.get_size());

//    for (int i = 0; i < this->array_size; ++i) {
//        this->p[i] = this->p[i] * a[i];
//    }
//}

//void MyArray::divideByElementwise(MyArray a) {
//    assert(this->array_size == a.get_size());

//    for (int i = 0; i < this->array_size; ++i) {
//        this->p[i] = this->p[i] / a[i];
//    }
//}

//void MyArray::multiplyElementwise(MyArray a, MyArray &res) {
//    assert(this->array_size == a.get_size());

//    for (int i = 0; i < this->array_size; ++i) {
//        res[i] = this->p[i] * a[i];
//    }
//}

//void MyArray::divideByElementwise(MyArray a, MyArray &res) {
//    assert(this->array_size == a.get_size());

//    for (int i = 0; i < this->array_size; ++i) {
//        res[i] = this->p[i] / a[i];
//    }
//}

bool MyArray::equalsToArray(MyArray a, float eps) {
    if (this->array_size != a.get_size()) {
        std::cout << "Dimension mismatch!\n";
        return false;
    }

    for (int i = 0; i < this->array_size; ++i) {
        if (std::abs(this->p[i] - a[i]) >= eps) {
            std::cout << "Array element " << i << " mismatch!\n";
            std::cout << this->p[i] << " vs " << a[i] << std::endl;
            return false;
        }
    }

    return true;
}

void MyArray::Hessian(float (*f)(MyArray x), Matrix &matr, float eps) {
	float tau = 0.1 * sqrt(eps);
	int k = 0, n = array_size;
	MyArray tmp1(n), tmp2(n), tmp3(n), tmp4(n);
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i == j) {
				for (k = 0; k < n; k++) {
					if (k == j) {
						tmp1[k] = p[k] + tau;
						tmp2[k] = p[k] - tau;
					} else {
						tmp1[k] = p[k];
						tmp2[k] = p[k];
					}				
				}			
				matr(i, j) = (f(tmp1) - 2 * f(*this) + f(tmp2)) / (tau * tau);
			} else {
				for (k = 0; k < n; k++) {
					if (k == i) {
						tmp1[k] = p[k] + tau;
						tmp2[k] = p[k] + tau;
						tmp3[k] = p[k] - tau;
						tmp4[k] = p[k] - tau;	
					}
					if (k == j) {
						tmp1[k] = p[k] + tau;
						tmp2[k] = p[k] - tau;
						tmp3[k] = p[k] + tau;
						tmp4[k] = p[k] - tau;
					} else {
						tmp1[k] = p[k];
						tmp2[k] = p[k];
						tmp3[k] = p[k];
						tmp4[k] = p[k];
					}					
				}
				matr(i, j) = (f(tmp1) - f(tmp2) - f(tmp3) + f(tmp4)) / (4 * tau * tau);
				matr(j, i) = (f(tmp1) - f(tmp2) - f(tmp3) + f(tmp4)) / (4 * tau * tau);
			}	
		}
	}
}





void MyArray::WriteToFile() {
	fstream out;
	out.open("output.txt", fstream::out);

	if (array_size > 0) {
		for (int i = 0; i < array_size; i++) {
			out << p[i] << " ";
		}
		out << endl;	
	} else {
		cout << "Incorrect data! : WriteToFile" << endl;
	}
}

void MyArray::ReadFromFile() {
	fstream in;
	in.open("input.txt", fstream::in);
	
	if (array_size > 0) {
		for (int i = 0; i < array_size; i++) {
			in >> p[i];
		}	
	} else {
		cout << "Incorrect data! : ReadFromFile" << endl;
	}
}


///////////////////
//////MATRIX///////
///////////////////

Matrix::Matrix() {
	this->row = 0;
	this->col = 0;
    this->m = nullptr;
}

Matrix::Matrix(int n) {
	this->row = n;
	this->col = n;

	this->m = new float[row * col];

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {			
			m[j + i * col] = 0.0;
		}		
	}
}

Matrix::Matrix(int row, int col) {
	this->row = row;
	this->col = col;

	this->m = new float[row * col];

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {			
			m[j + i * col] = 0.0;
		}		
	}
}
	
Matrix::Matrix(const Matrix &a) {
	this->row = a.row;
	this->col = a.col;
		
	m = new float[row * col];

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {			
			m[j + i * col] = a.m[j + i * col];
		}		
	}
}
	

Matrix::~Matrix() {
	delete [] m;
}

void Matrix::Resize(int row, int col) {
    this->row = row;
    this->col = col;
    this->m = new float[row * col];

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            m[j + i * col] = 0.0;
        }
    }
}

float & Matrix::operator ()(int i, int j) {
	return m[j + i * col];
}

void Matrix::Show() {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {			
			cout << "[" << j + i * col << "] = " << m[j + i * col] << " ";
		}
		cout << endl;		
	}
	cout << endl;
}

int Matrix::CountNonzero() {
    int count = 0;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (!(m[j + i * col] == 0.0f)) ++count;
        }
    }

    return count;
}

bool Matrix::equalsToMatrix(Matrix &a, float eps) {
    if (col != a.get_col() || row != a.get_row()) {
        std::cout << "Dimension mismatch!\n";
        return false;
    }
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (std::abs(m[j + i * col] - a(j,i)) >= eps) {
                std::cout << "Matrix element (" << j << ", " << i << ") mismatch!\n";
                std::cout << m[j + i * col] << " vs " << a(j,i) << std::endl;
                return false;
            }
        }
    }

    return true;
}

void Matrix::Set(int index1, int index2, float value) {
	if ((index1 >= 0) && (index2 >= 0) && (index1 < row) && (index2 < col)) {
		m[index2 + index1 * col] = value;
    } else {
        fprintf(stderr,"Indexes exceed matrix dimensions!\n");
        exit(1);
    }
}

float Matrix::Get(int index1,int index2)
{
    if ((index1 >= 0) && (index2 >= 0) && (index1 < row) && (index2 < col)) {
        return m[index2 + index1 * col];
    } else {
        fprintf(stderr,"Indexes exceed matrix dimensions!\n");
        exit(1);
    }
}

int Matrix::get_row() {
	return row;
}

int Matrix::get_col() {
	return col;
}

float* Matrix::get_data() {
    return m;
}

void Matrix::zap() {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m[j + i * col] = rand() % 10;
		}
	}
}
	
Matrix Matrix::operator =(const Matrix &a) {
	this->row = a.row;
	this->col = a.col;

	this->m = new float[row * col];

    for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
      			m[j + i * col] = a.m[j + i * col];
		}
	}
	return *this;
}

void Matrix::scale(float value) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m[j + i * col] *= value;
		}
	}
}

Matrix Matrix::Sum(Matrix &a) {
    Matrix temp(a.row, a.col);
    if (row == a.row || col == a.col) {
        for (int i = 0; i < temp.row; i++) {
            for (int j = 0; j < temp.col; j++) {
                temp(i, j) = m[j + i * col] + a(i, j);
            }
        }
    } else {
        cout << "ERROR! Sum: row1!=row2 or col1!=col2" << endl;
    }
    return temp;
}

Matrix Matrix::weightedSum(Matrix &a, float alpha, float beta) {
    Matrix temp(a.row, a.col);
    if (row == a.row || col == a.col) {
        for (int i = 0; i < temp.row; i++) {
            for (int j = 0; j < temp.col; j++) {
                temp(i, j) = alpha * m[j + i * col] + beta * a(i, j);
            }
        }
    } else {
        cout << "ERROR! Sum: row1!=row2 or col1!=col2" << endl;
    }
    return temp;
}

Matrix Matrix::Difference(Matrix &a) {
	Matrix temp(a.row, a.col);
    if (row == a.row || col == a.col) {
		for (int i = 0; i < temp.row; i++) {
			for (int j = 0; j < temp.col; j++) {			
				temp(i, j) = m[j + i * col] - a(i, j);
			}		
		}
	} else {
		cout << "ERROR! Difference: row1!=row2 or col1!=col2" << endl;
	}
	return temp;
}
	
Matrix Matrix::Product(Matrix &a) {
	Matrix temp(row, a.col);
	if (col == a.row) {	
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < a.col; j++) {
				temp(i, j) = 0.0;
				for (int k = 0; k < col; k++) {
					temp(i, j) += m[k + i * col] * a(k, j);
				}				
			}
		}
	} else {
        cout << "ERROR! Product: col1 != row2" << endl;
	}
	return temp;
}

MyArray Matrix::Product(MyArray &v) {
    MyArray temp(row);
    if (col == v.get_size()) {
        for (int i = 0; i < row; i++) {
            temp[i] = 0.0;
            for (int j = 0; j < col; j++) {
                temp[i] += m[j + i * col] * v[j];
            }
        }
    } else {
       cout << "ERROR! Product: col != vector size" << endl;
    }
    return temp;
}
	
void Matrix::LU_decomposition(Matrix &L, Matrix &U, int n) {
	U = *this;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			L.m[i + j * col] = U.m[i + j * col] / U.m[i + i * col];
		}	
	}

	for (int k = 1; k < n; k++) {
		for (int i = k - 1; i < n; i++) {
			for (int j = i; j < n; j++) {
				L.m[i + j * col] = U.m[i + j * col] / U.m[i + i * col];
			}
		}
	
		for (int i = k; i < n; i++) {
			for (int j = k - 1; j < n; j++) {
				U.m[j + i * col] = U.m[j + i * col] - L.m[k - 1 + i * col] * U.m[j + (k - 1) * col];
			}
		}	
	}
	L.isDiag = true;
	U.isDiag = true;
}

Matrix & Matrix::transpose() {
	float temp;
	if (col == row) {		
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				if (i > j) {
					temp = m[j + i * col];
					m[j + i * col] = m[i + j * row];
					m[i + j * row] = temp;
				}
			}
		}
        this->isSym = true;
	} else {
		cout << "ERROR! transpose: col != row" << endl;
	}	
    return *this;
}

void Matrix::transpose2() {
	Matrix temp_matr(row, col);
	temp_matr = *this;

	int temp_size = row;
	this->row = col;
	this->col = temp_size;	

	if (col == row) {		
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				//temp_matr(i, j) = m[i + j * row];
				m[j + i * col] = temp_matr.m[i + j * row];
			}
		}
        //this->isSym = true;
	} else {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
                m[j + i * col] = temp_matr.m[i + j * row];
			}
		}
	}	
    //return *this;
}

void Matrix::transpose_not_square() {

}

//not ideal
Matrix & Matrix::Gauss() {
	float t = 0.0;
	for (int k = 0; k < row; k++) {
		for (int i = k + 1; i < row; i++) {
			t = (float)m[k + i * col] / m[k + k * col];
			for (int j = 0; j < col; j++) {
				m[j + i * col] -= m[j + k * col] * t;
			}
		}
	}
	this->isDiag = true;
	return *this;
}

//not ideal
float Matrix::det_gauss() {
	float determinant = 1.0;
	Matrix temp(row, col);
	temp = *this;
	
	if (!(temp.isDiag)) 
		temp.Gauss();
	
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (i == j) {
				determinant *= temp.m[j + i * col];
			}
		}
	}
	return determinant;
}

void Matrix::Get_matrix(int n, Matrix &temp_matr, int indRow, int indCol) {
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr(ki, kj) = m[j + i * col];
					kj++;
				}
			}
			ki++;
		}
	}
}

float det3x3(float a0, float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
    Matrix temp(3, 3);
    temp(0, 0) = a0;
    temp(0, 1) = a1;
    temp(0, 2) = a2;
    temp(1, 0) = a3;
    temp(1, 1) = a4;
    temp(1, 2) = a5;
    temp(2, 0) = a6;
    temp(2, 1) = a7;
    temp(2, 2) = a8;
    return temp.det(3);
}

float Matrix::det(int n) {
	float temp = 0;
	int k = 1;
	if (n == 1) {
		return m[0];
	} else if (n == 2) {
		return m[0 + 0 * col] * m[1 + 1 * col] - m[0 + 1 * col] * m[1 + 0 * col];
    } else if (n == 3) {
        return m[0]*m[4]*m[8] +
                m[1]*m[6]*m[5] +
                m[2]*m[3]*m[7] -
                m[6]*m[4]*m[2] -
                m[0]*m[5]*m[7] -
                m[1]*m[3]*m[8];
    } else if (n == 4) {
        float v1 = det3x3(m[5], m[6], m[7], m[9], m[10], m[11], m[13], m[14], m[15]);
        float v2 = det3x3(m[1], m[2], m[3], m[9], m[10], m[11], m[13], m[14], m[15]);
        float v3 = det3x3(m[1], m[2], m[3], m[5], m[6], m[7], m[13], m[14], m[15]);
        float v4 = det3x3(m[1], m[2], m[3], m[5], m[6], m[7], m[9], m[10], m[11]);
        return v1 - v2 + v3 - v4;
    } else if (n >= 5) {
		for (int i = 0; i < n; i++) {
			int p = n - 1;
			Matrix temp_matr(p, p);
			this->Get_matrix(n, temp_matr, 0, i);
			temp += k * m[i + 0 * col] * temp_matr.det(p);
            k -= k;
		}
		return temp;
	}
    cout << "Det = 0.0!!!";
    return 0.0;
}

void Matrix::inverse(Matrix &matr, int n, bool isDiag) {
	Matrix inverse_matr(n, n);
    float determinant = this->det(n);
	


	if (determinant) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int p = n - 1;
				Matrix temp_matr(p, p);
				this->Get_matrix(n, temp_matr, i, j);
                inverse_matr(i, j) = pow(-1.0, i + j + 2) * temp_matr.det(p) / determinant;
			}
		}
	} else {
		cout << "ERROR! inverse: DETERMINANT = 0!" << endl;
		exit(-1);
	}
	inverse_matr.transpose();
	matr = inverse_matr;
}

void Matrix::Solve_Gauss_reverse(Matrix matr, MyArray B, MyArray &res, int n, bool isDown) {
	if (!isDown) {
		for (int i = n - 1; i >= 0; i--) {
			float temp = 0.0;
			for (int j = i + 1; j < n; j++) {
				temp += matr(i, j) * res[j];
			}
		res[i] = (B[i] - temp) / matr(i, i);
		}
	} else {
		for (int i = 0; i < n; i++) {
			float temp = 0.0;
			for (int j = i - 1; j >= 0; j--) {
				temp += matr(i, j) * res[j];
			}
		res[i] = (B[i] - temp) / matr(i, i);
		}
	}
} 


void Matrix::LU_solve(MyArray B, MyArray &result, int n) {
	MyArray res_tmp(n);
	Matrix L(n, n), U(n, n);
	
    this->LU_decomposition(L, U, n);
    std::cout << "LU_decomposition success\n";
	
	// std::cout << "L = \n";
	// L.Show();
	// std::cout << "U = \n";
	// U.Show();

	Solve_Gauss_reverse(L, B, res_tmp, n, 1);
	Solve_Gauss_reverse(U, res_tmp, result, n, 0);
}

void Matrix::CGM_solve(Matrix A, MyArray B, MyArray &x_k, int n) {
	int k = 1;

    float eps = 0.00001;
	float *z_k = new float[n];
	float *r_k = new float[n];
  	float *Az = new float[n]; 
	float alpha, beta, mf = 0.0;
	float Spr, Spr1, Spz;

  for (int i = 0; i < n; i++) {
    mf += B[i] * B[i];
    x_k[i] = 0.2;
  }

  for (int i = 0; i < n; i++) {
    Az[i] = 0.0;
    for (int j = 0; j < n; j++) {
      Az[i] += A(i, j) * x_k[j];
    }
    r_k[i] = B[i] - Az[i];
    z_k[i] = r_k[i];
  }

  do{
    Spz=0.0;
    Spr=0.0;
    for (int i = 0; i < n; i++) {
      Az[i] = 0.0;
      for (int j = 0; j < n; j++) {
        Az[i] += A(i, j) * z_k[j];   
      }
      Spz += Az[i] * z_k[i];
      Spr += r_k[i] * r_k[i];
    }
    alpha = Spr / Spz;

    Spr1 = 0.0;
    for (int i = 0; i < n; i++) {
      x_k[i] += alpha * z_k[i];
      r_k[i] -= alpha * Az[i];
      Spr1 += r_k[i] * r_k[i];
      cout << "Iter #" << k;
      cout << " " << "X[" << i << "] = " << x_k[i] << endl;
    }
    cout << endl;
    k++;

    beta = Spr1 / Spr;
  
    for (int i = 0; i < n; i++) {
      z_k[i] = r_k[i] + beta * z_k[i];
    }
  } while(Spr1 / mf > eps * eps);

  cout << endl;

  delete [] Az;
  delete [] z_k;
  delete [] r_k;
}


void Matrix::WriteToFile() {
	fstream out;
	out.open("output.txt", fstream::out);

	if (row > 0 && col > 0) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				out << m[j + i * col] << " ";
			}
			out << endl;
		}
	} else {
		cout << "Incorrect data! : WriteToFile" << endl;
	}
}

void Matrix::ReadFromFile() {
	fstream in;
	in.open("input.txt", fstream::in);
	if (row > 0 && col > 0) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				in >> m[j + i * col];
			}
		}
	} else {
		cout << "Incorrect data! : ReadFromFile" << endl;
	}
}


///////////////////////
/////SPARSE MATRIX/////
///////////////////////

SparseMatrixCOO::SparseMatrixCOO() {
    this->sparse_size = 0;
    this->x = nullptr;
    this->y = nullptr;
    this->data = nullptr;
    this->sorted = false;
    this->data_pointer = 0;
}

SparseMatrixCOO::SparseMatrixCOO(int size) {
    this->sparse_size = size;
    this->x = size == 0? nullptr : new int[sparse_size];
    this->y = size == 0? nullptr : new int[sparse_size];
    this->data = size == 0? nullptr : new float[sparse_size];
    this->sorted = false;
    this->data_pointer = 0;
}

SparseMatrixCOO::SparseMatrixCOO(const SparseMatrixCOO &m) {
    sparse_size = m.sparse_size;
    this->x = sparse_size == 0? nullptr : new int[sparse_size];
    this->y = sparse_size == 0? nullptr : new int[sparse_size];
    this->data = sparse_size == 0? nullptr : new float[sparse_size];
    this->sorted = false;
    this->data_pointer = sparse_size;

    for (int i = 0; i < sparse_size; i++) {
        x[i] = m.x[i];
        y[i] = m.y[i];
        data[i] = m.data[i];
    }
}

SparseMatrixCOO SparseMatrixCOO::operator =(const SparseMatrixCOO &m) {
    sparse_size = m.sparse_size;
    this->x = new int[sparse_size];
    this->y = new int[sparse_size];
    this->data = new float[sparse_size];
    this->data_pointer = sparse_size;

    for (int i = 0; i < sparse_size; i++) {
        x[i] = m.x[i];
        y[i] = m.y[i];
        data[i] = m.data[i];
    }
    return *this;
}

SparseMatrixCOO::~SparseMatrixCOO() {
	delete [] x;
	delete [] y;
	delete [] data;
}

void SparseMatrixCOO::resize() {
    CheckRunTime(__func__)
    int *x_new = new int[sparse_size];
    int *y_new = new int[sparse_size];
    float *data_new = new float[sparse_size];
    for (int i = 0; i < sparse_size; i++) {
        x_new[i] = x[i];
        y_new[i] = y[i];
        data_new[i] = data[i];
    }
    delete [] x;
    delete [] y;
    delete [] data;
    x = x_new;
    y = y_new;
    data = data_new;
}

void SparseMatrixCOO::resize(int size) {
    assert(size > sparse_size);
    int *x_new = new int[size];
    int *y_new = new int[size];
    float *data_new = new float[size];
    for (int i = 0; i < size; i++) {
        x_new[i] = x[i];
        y_new[i] = y[i];
        data_new[i] = data[i];
    }
    delete [] x;
    delete [] y;
    delete [] data;
    x = x_new;
    y = y_new;
    data = data_new;
    sparse_size = size;
}

int SparseMatrixCOO::get_size() { return sparse_size; }

void SparseMatrixCOO::set_value(int row, int col, float value) {
//    for (int i = 0; i < sparse_size; i++) {
//        if (x[i] == row && y[i] == col) {
//            data[i] = value;
//        }
//    }
    for (int i = 0; i < sparse_size; ++i) {
        if (x[i] == row && y[i] == col) {
            data[i] = value;
        }
    }

}

void SparseMatrixCOO::write_value(int row, int col, float value) {
    assert(this->data_pointer < this->sparse_size);

    if (!(value == 0.0f)) {
        this->x[this->data_pointer] = row;
        this->y[this->data_pointer] = col;
        this->data[this->data_pointer] = value;

        ++data_pointer;
    }

}

void SparseMatrixCOO::pop(int row, int col)
{
    for (int i = 0; i < sparse_size; i++) {
        if (x[i] == row && y[i] == col) {
            x[i] = x[sparse_size-1];
            y[i] = y[sparse_size-1];
            data[i] = data[sparse_size-1];
            --sparse_size;
            break;
        }
    }
}

int SparseMatrixCOO::get_x(int index) {
    return x[index];
}

int SparseMatrixCOO::get_y(int index) {
    return y[index];
}

float SparseMatrixCOO::get_value(int index) {
    return data[index];
}

int* SparseMatrixCOO::get_x() {
    return x;
}

int* SparseMatrixCOO::get_y() {
    return y;
}

float* SparseMatrixCOO::get_data() {
    return data;
}

SparseMatrixCOO SparseMatrixCOO::DeleteZeros() {
    CheckRunTime(__func__)
    int new_size = sparse_size;
    int k = 0;
    SparseMatrixCOO temp(sparse_size);
    for (int i = 0; i < sparse_size; ++i) {
        if (data[i] != 0) {
            temp.x[k] = x[i];
            temp.y[k] = y[i];
            temp.data[k] = data[i];
            ++k;
        } else {
            --new_size;
        }
    }
    cout<<"new_size = " << new_size<<endl;
    temp.resize(new_size);
    return temp;
}

void SparseMatrixCOO::DeleteZerosOrder() {
    CheckRunTime(__func__)
    int new_size = sparse_size;
    for (int i = 0, k = 0; i < sparse_size; ++i) {
        if (data[i] != 0.0) {
            x[k] = x[i];
            y[k] = y[i];
            data[k] = data[i];
            ++k;
        } else {
            --new_size;
        }
    }
    std::cout << "new_size = " << new_size << std::endl;
    resize(new_size);
}

int SparseMatrixCOO::CountNonZero() {
    int nonzero = 0;
    for (int i = 0; i < sparse_size; i++) {
        if (data[i] != 0) {
            nonzero++;
        }
    }
    return nonzero;
}

//int CountNonZero(std::vector<Triplet> t) {
//    int nonzero = 0;
//    for (int i = 0; i < t.size(); i++) {
//        if (t[i].get_value() != 0.0) {
//            nonzero++;
//        } else {
//            t.erase(t.begin() + i);
//        }
//    }
//    return nonzero;
//}


void SparseMatrixCOO::ConvertTripletToSparse(std::vector<Triplet> t) {
    CheckRunTime(__func__)
    int new_size = sparse_size;
    int k = 1;
    int index = 0;
    bool changeme = false;

    x[0] = t[0].x_value;
    y[0] = t[0].y_value;
    data[0] = t[0].value;

    for (int i = 1; i < sparse_size; i++ ) {
        for (int j = 0; j < k; j++) {
            if (t[i].x_value == x[j] && t[i].y_value == y[j]) {
                changeme = true;
                index = j;
                break;
            }
        }
        //index = FindIndex(k, t[i].x_value, t[i].y_value, x, y, sparse_size);

        if (changeme == true) {
            data[index] += t[i].value;
            new_size--;
            changeme = false;
        } else {
            x[k] = t[i].x_value;
            y[k] = t[i].y_value;
            data[k] = t[i].value;
            k++;
        }
   }
   sparse_size = new_size;



//    bool changeme = false;
//    for (int i = 0; i < sparse_size; i++ ) {
//        if (i == 0) {
//            x[i] = t[i].x_value;
//            y[i] = t[i].y_value;
//            data[i] = t[i].value;
//            this->v.push_back(t[i]);
//        } else {
//            for (int j = 0; j < i; j++) {
//                if (t[i].x_value == x[j] && t[i].y_value == y[j]) {
//                    changeme = true;
//                    t[i].Show();
//                }
//            }
//            if (changeme == true) {
//                for (int j = 0; j < i; j++) {
//                    if (t[i].x_value == x[j] && t[i].y_value == y[j]) {
//                        data[j] += t[i].value;
//                        new_size--;
//                        changeme = false;
//                    }
//                }
//            } else {
//                x[i] = t[i].x_value;
//                y[i] = t[i].y_value;
//                data[i] = t[i].value;
//                this->v.push_back(t[i]);
//            }
//        }
//    }
//    sparse_size = new_size;

}

// void SparseMatrixCOO::ConvertToMatrix(Matrix &M) {
// 	int k = 0;
// 	for (int i = 0; i < M.get_row(); i++) {
// 		for (int j = 0; j < M.get_col(); j++) {
// 			if (x[k] == i && y[k] == j) {
// 				M(i, j) = data[k];
// 				//cout<<data[k]<<" ";
// 				k++;
// 			} else {
// 				M(i, j) = 0.0;
// 			}
// 		}
// 	}
// }

void SparseMatrixCOO::ConvertToMatrix(Matrix &M) {
    // Make sure you set the right dimensions for matrix M!!!
    int n = M.get_row();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
			for (int k = 0; k < sparse_size; k++) {
				if (i == x[k] && j == y[k] && data[k] != 0) {
					M(i, j) += data[k];	
				}
			}
		}
	}
}

MyArray SparseMatrixCOO::MultiplyByVector(MyArray v, int res_size) {
    MyArray res(res_size);
    for (int i = 0; i < this->sparse_size; ++i) {
        res[this->x[i]] += this->data[i] * v[this->y[i]];
    }

    return res;
}

MyArray SparseMatrixCOO::MultiplyTransposedByVector(MyArray v, int res_size) {
    MyArray res(res_size);
    for (int i = 0; i < this->sparse_size; ++i) {
        res[this->y[i]] += this->data[i] * v[this->x[i]];
    }

    return res;
}

//void SparseMatrixCOO::MultiplyByVector(MyArray v, MyArray &res) {
//    for (int i = 0; i < this->sparse_size; ++i) {
//        res[this->x[i]] += this->data[i] * v[this->y[i]];
//    }
//}

//void SparseMatrixCOO::MultiplyTransposedByVector(MyArray v, MyArray &res) {
//    for (int i = 0; i < this->sparse_size; ++i) {
//        res[this->y[i]] += this->data[i] * v[this->x[i]];
//    }
//}

void SparseMatrixCOO::SortIt() {
    CheckRunTime(__func__)
	int temp;
	float temp_value;
    for (int i = 0; i < sparse_size; i++) {
        for (int j = 0; j < sparse_size - 1; j++) {
            if (x[j] > x[j + 1]) {
				temp = x[j];
				x[j] = x[j + 1];
				x[j + 1] = temp;

				temp = y[j];
				y[j] = y[j + 1];
				y[j + 1] = temp;

				temp_value = data[j];
				data[j] = data[j + 1];
				data[j + 1] = temp_value;
			}
		}
	}

    for (int i = 0; i < sparse_size; i++) {
        for (int j = 0; j < sparse_size - 1; j++) {
            if (y[j] > y[j + 1] && x[j] == x[j + 1]) {
                temp = x[j];
                x[j] = x[j + 1];
                x[j + 1] = temp;

                temp = y[j];
                y[j] = y[j + 1];
                y[j + 1] = temp;

                temp_value = data[j];
                data[j] = data[j + 1];
                data[j + 1] = temp_value;
            }
        }
    }
    this->sorted = true;
}

bool SparseMatrixCOO::isSorted() {
    return this->sorted;
}

void SparseMatrixCOO::ConvertToCSR(int *ptr, int *ind, float *data_csr, int n) {
        for (int i = 0; i < n + 1; i++)
        {
            //ptr[i] = 0.0;
        }

        for (int i = 0; i < sparse_size; i++)
        {
            data_csr[i] = data[i];
            ind[i] = y[i];
            //ptr[x[i] + 1]++;
        }

        for (int i = 0; i < n; i++)
        {
            //ptr[i + 1] += ptr[i];
        }
}

void SparseMatrixCOO::ShowData() {
    for (int i = 0; i < sparse_size; i++) {
        cout << data[i] << " ";
    }
}

void SparseMatrixCOO::SparseLU() {
//    for (int i = 0; i < n; i++) {
//        for (int j = i; j < n; j++) {
//            L.m[i + j * col] = U.m[i + j * col] / U.m[i + i * col];
//        }
//    }

//    for (int k = 1; k < n; k++) {
//        for (int i = k - 1; i < n; i++) {
//            for (int j = i; j < n; j++) {
//                L.m[i + j * col] = U.m[i + j * col] / U.m[i + i * col];
//            }
//        }

//        for (int i = k; i < n; i++) {
//            for (int j = k - 1; j < n; j++) {
//                U.m[j + i * col] = U.m[j + i * col] - L.m[k - 1 + i * col] * U.m[j + (k - 1) * col];
//            }
//        }
//    }
//        float value1;
//        int x1, y1;
//    cout << "SparseLU ------\n";
//    SparseMatrixCOO U(*this);
//    SparseMatrixCOO L(sparse_size);
//        for (int i = 0; i < sparse_size; i++) {
//            if (x[i] == y[i]) {
//                value1 = data[i];
//                x1 = x[i];
//                y1 = y[i];
//            }
//            for (int j = i; j < sparse_size; j++) {
//                L.data[j] = U.data[j] / value1;
//            }
//        }
//            for (int k = 1; k < sparse_size; k++) {
//                for (int i = k - 1; i < sparse_size; i++) {
//                    if (x[i] == y[i]) {
//                        value1 = data[i];
//                        x1 = x[i];
//                        y1 = y[i];
//                    }
//                    for (int j = i; j < sparse_size; j++) {
//                        L.data[j] = U.data[j] / value1;
//                    }
//                }

//                for (int i = k; i < sparse_size; i++) {
//                    for (int j = k - 1; j < sparse_size; j++) {
//                        U.data[j] = U.data[j] - L.data[k - 1 + i * col] * U.data[j + (k - 1) * col];
//                    }
//                }
//            }
//    cout<<"U=\n";
//    U.Show();
//    cout<<"L=\n";
//    L.Show();



//    for (int i = 0; i < sparse_size; i++) {
//        if (x[i] == y[i]) {
//            value1 = data[i];
//            x1 = x[i];
//            y1 = y[i];
//        }
//        for (int j = 0; j < sparse_size; j++) {
//             if (x[j] == x1) {
//                data[j] = (data[j] * value1 > 0) ? data[j] - data[j] / value1 * data[j]
//                                                 : data[j] + data[j] / value1 * data[j];
//            } else if (y[j] == y1) {

//            }

//        }
//    }
}

void SparseMatrixCOO::set_diag_elements() {
    if (!this->sorted) {
        this->SortIt();
    }
    for (int i = 0; i < sparse_size; ++i) {
        if (x[i] == y[i]) {
            this->diag_elems.push_back(data[i]);
        }
    }
}

std::vector<float> SparseMatrixCOO::get_diag_elements() {
    assert(diag_elems.size()!=0);
    return this->diag_elems;
}

void SparseMatrixCOO::CGM_solve(MyArray B, MyArray &x_k, float eps) {
    CheckRunTime(__func__)
    int n = B.get_size();
    int k = 1;

    float *z_k = new float[n];
    float *r_k = new float[n];
    float *Az = new float[n];
    float alpha, beta, mf = 0.0;
    float Spr, Spr1, Spz;

  for (int i = 0; i < n; i++) {
    mf += B[i] * B[i];
    x_k[i] = 0.0;
  }

  for (int i = 0; i < n; i++) {
    Az[i] = 0.0;
  }
  for (int i = 0; i < sparse_size; i++) {
      Az[x[i]] += data[i] * x_k[y[i]];
  }
  for (int i = 0; i < n; i++) {
    r_k[i] = B[i] - Az[i];
    z_k[i] = r_k[i];
  }

  do{
    Spz=0.0;
    Spr=0.0;
    for (int i = 0; i < n; i++) {
      Az[i] = 0.0;
    }
    for (int i = 0; i < sparse_size; i++) {
        Az[x[i]] += data[i] * z_k[y[i]];
    }
    for (int i = 0; i < n; ++i) {
      Spz += Az[i] * z_k[i];
      Spr += r_k[i] * r_k[i];
    }
    alpha = Spr / Spz;

    Spr1 = 0.0;
    for (int i = 0; i < n; ++i) {
      x_k[i] += alpha * z_k[i];
      r_k[i] -= alpha * Az[i];
      Spr1 += r_k[i] * r_k[i];
      //cout << "Iter #" << k;
      //cout << " " << "X[" << i << "] = " << x_k[i] << endl;
    }
    //cout << endl;
    k++;

    beta = Spr1 / Spr;

    for (int i = 0; i < n; ++i) {
      z_k[i] = r_k[i] + beta * z_k[i];
    }
  } while(Spr1 / mf > eps * eps);

  //cout << endl;

  delete [] Az;
  delete [] z_k;
  delete [] r_k;
}

void SparseMatrixCOO::PCG_solve(MyArray B, MyArray &x_k, float eps) {
    CheckRunTime(__func__)
    int n = B.get_size();
    int k = 1;

    float *z_k = new float[n];
    float *r_k = new float[n];
    float *Az = new float[n];
    float alpha, beta, mf = 0.0;
    float Spr, Spr1, Spz, Spr_initial = 0.0;

  for (int i = 0; i < n; i++) {
    mf += B[i] * B[i]; // change stop condition
    x_k[i] = 0.0;
  }

  for (int i = 0; i < n; i++) {
    Az[i] = 0.0;
  }
  for (int i = 0; i < sparse_size; i++) {
      Az[x[i]] += data[i] * x_k[y[i]];
  }
  for (int i = 0; i < n; i++) {
    r_k[i] = B[i] - Az[i];
    z_k[i] = r_k[i] / this->diag_elems[i];
    Spr_initial += r_k[i] * r_k[i] / this->diag_elems[i];
  }

  std::cout.precision(16);
  std::cout << "gamma0\t\t= " << Spr_initial << std::endl;

  int n_iter = 0;
  do{
    ++n_iter;
    std::cout << "Iteration #" << n_iter << std::endl;
    Spz=0.0;
    Spr=0.0;
    for (int i = 0; i < n; i++) {
      Az[i] = 0.0;
    }
    for (int i = 0; i < sparse_size; i++) {
        Az[x[i]] += data[i] * z_k[y[i]];
    }
    for (int i = 0; i < n; ++i) {
      Spz += Az[i] * z_k[i];
      Spr += r_k[i] * r_k[i] / this->diag_elems[i];
    }
    alpha = Spr / Spz;
    std::cout << "sumElem\t\t\t= " << Spz << std::endl;
    std::cout << "alpha\t\t\t= " << alpha << std::endl;
    std::cout << "alpha numerator (gamma)\t= " << Spr << std::endl;
    std::cout << "alpha denominator\t= " << Spz << std::endl;

    Spr1 = 0.0;
    for (int i = 0; i < n; ++i) {
      x_k[i] += alpha * z_k[i];
      r_k[i] -= alpha * Az[i];
      Spr1 += r_k[i] * r_k[i] / this->diag_elems[i];
      //cout << "Iter #" << k;
      //cout << " " << "X[" << i << "] = " << x_k[i] << endl;
    }
    //cout << endl;
    k++;

    beta = Spr1 / Spr;

    std::cout << "beta\t\t\t= " << beta << std::endl;
    std::cout << "gamma_new\t\t= " << Spr1 << std::endl;

    for (int i = 0; i < n; ++i) {
      z_k[i] = r_k[i] / this->diag_elems[i] + beta * z_k[i];
    }
    std::cout << std::endl;
//  } while(Spr1 / mf > eps * eps);
    } while(Spr1 > eps * Spr_initial);

  //cout << endl;

  delete [] Az;
  delete [] z_k;
  delete [] r_k;
}


void SparseMatrixCOO::Show() {
	cout << "X\tY\tValue\n";
    for (int i = 0; i < sparse_size; ++i) {
		cout << x[i] << "\t" << y[i] << "\t" << data[i] << endl;
	}
}

// Warning: the ShowAsMatrix* methdods work well only if the sparse matrix is sorted (call method .SortIt())
void SparseMatrixCOO::ShowAsMatrixNumber(int start, int end, int n) {
    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            bool cond_start_end = (i >= start && i <= end && j >= start && j <= end);
            if (x[k] == i && y[k] == j) {
                if (cond_start_end)
                    cout << data[k] << " ";
                ++k;
            } else {
                if (cond_start_end)
                    cout << "- ";
            }
        }
        cout << endl;
    }
}

// Warning: see warning before
void SparseMatrixCOO::ShowAsMatrixSymbol(int start, int end, int n) {
    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            bool cond_start_end = (i >= start && i <= end && j >= start && j <= end);
            if (x[k] == i && y[k] == j) {
                if (cond_start_end)
                    cout << "A ";
                ++k;
            } else {
                if (cond_start_end)
                    cout << "- ";
            }
        }
        cout << endl;
    }
}

// Warning: see warning before
void SparseMatrixCOO::ShowAsMatrix(int start, int end, int n) {
    ShowAsMatrixSymbol(start, end, n);
}
