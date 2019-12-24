#include <iostream>
#include <cmath>

using namespace std;

const double eps = 0.00000001;

class Matrix {
    private:
        int rows;
        int cols;
        int alloc_rows;
        float** matrix;
        bool expanded;
    public:
        Matrix(int r, int c) {
            rows = r;
            cols = c;
            matrix = new float* [rows];
            alloc_rows = r;
            for (int i = 0; i < rows; ++i) {
                matrix[i] = new float[cols];
                if(matrix[i] == nullptr) {
                    deleteMatrix(i);
                }
            }
        }
        ~Matrix() {
            deleteMatrix(alloc_rows);
        }
        void swap_rows(int row1, int row2) {
            for(int i=0; i < cols; i++) {
                float tmp = matrix[row1][i];
                matrix[row1][i] = matrix[row2][i];
                matrix[row2][i] = tmp;
            }
        }
        int set_expansion(bool exp) {
            expanded = exp;
        }
        int get_rows(){
            return rows;
        }
        float* get_row(int i){
            return matrix[i];
        }
        float el(int row, int col){
            return matrix[row][col];
        }
        int working_cols() {
            return expanded ? cols - 1 : cols;
        }
        int get_cols() {
            return cols;
        }
        void deleteMatrix(int rows) {
            for(int i=0;i < rows; i++) {
                delete [] matrix[i];
            }
            delete [] matrix;
        }
        bool rowIsEmpty(int row) {
            for(int i = 0; i < cols;i++) {
                if(fabs(matrix[row][i]) > eps) {
                    return false;
                }
            }
            return true;
        }
        int removeRow(int row) {
            if(row == rows - 1) {
                delete[] matrix[row];
                rows--;
                return rows;
            }

            for(int i = row; i < rows - 1; i++) {
                for(int j = 0; j < cols; j++) {
                    matrix[i][j] = matrix[i + 1][j];
                }
            }
            delete [] matrix[rows - 1];
            rows--;
            return rows;
        }
        int devideRow(int row, float devider) {
            if(fabs(devider - 1) > eps) {
                for(int i = 0; i < cols; i++) {
                    if(matrix[row][i] != 0) {
                        matrix[row][i] /= devider;
                    }   
                }
                return 1;
            }
            return -1;
        }
        void printMatrix() {
            cout << endl;
            for(int i=0; i < rows;i++) {
                for(int j=0;j < cols;j++) {
                    cout << matrix[i][j] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
        void printEquations() {
            cout << endl;
            if(!expanded) {
                return;
            }
            
            for(int i = 0; i < rows; i++) {
                int working = working_cols();
                for(int j = 0; j < working ; j++) {
                    if(matrix[i][j] != 0) {
                        cout << (matrix[i][j] > 0 ? " + " : " - ") << fabs(matrix[i][j]) << "*x" << j + 1;
                    }
                }
                cout << " = " << matrix[i][cols - 1] << endl;
            }
        }
        void transformRows(int heading_elem_row, int heading_elem_col) {
            for(int i = 0; i < rows; i++) {
                if(i == heading_elem_row) continue;
                float mult = matrix[i][heading_elem_col] / matrix[heading_elem_row][heading_elem_col];
                if(mult != 0) {
                    cout << "R" << i << (mult > 0 ? " - " : " + " ) << fabs(mult) << "R" << heading_elem_row << endl;
                    for(int j = 0; j < cols; j++) {
                        matrix[i][j] -= matrix[heading_elem_row][j] * mult;
                    }
                }
            }
        }
        void readMatrix(){
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++ ) {
                    cin >> matrix[i][j];
                }
            }
        }
};


class Method {
    protected:
        Matrix* matrix;
        float gcd(float a, float b) { 
            if (a < b) {
                return gcd(fabs(b), fabs(a));
            }
            if (fabs(b) < eps) {
                return a; 
            } else {
                return (gcd(fabs(b), fabs(a - floor(a / b) * b)));         
            }
        } 

        float findGCD(int row) { 
            int colums = matrix -> get_cols(); 
            float* arr = matrix -> get_row(row);
            float result = arr[0];

            for (int i = 1; i < colums; i++) 
            { 
                result = gcd(fabs(arr[i]), fabs(result)); 
        
                if(fabs(result) - 1 < eps) 
                { 
                return 1; 
                } 
            } 
            return result; 
        }
        void removeEmptyRows() {
            int rows = matrix -> get_rows();
            for(int i = 0; i < rows; i++) {
                if(matrix -> rowIsEmpty(i)) {
                    matrix -> removeRow(i);
                }
            }
        }
        void tryToDevideRows() {
            int rows = matrix -> get_rows();
            for(int i = 0; i < rows; i++) {
                float devider = findGCD(i);

                if(fabs(devider - 1) > eps) {
                    matrix -> devideRow(i, devider);
                }
            }
        }
    public:
        Method(Matrix* m) {
            matrix = m;
        }
        virtual void setHeadingElem() = 0;
        virtual void solution() = 0; 
};
class GaussJordan: public Method {
    private:
        float heading_elem;
        int heading_col;
        int heading_row;
    public:    
        GaussJordan(Matrix* m) : Method(m) {
            heading_elem = 0;
            heading_col = -1;
            heading_row = -1;
        };

        void setHeadingElem() {
            heading_row += 1;
            heading_col += 1;
            heading_elem = matrix -> el(heading_row, heading_col);
        }
        int searhForSwapRow() {
            int rows = matrix -> get_rows();
            for(int i = heading_row + 1; i < rows; i++) {
                if(matrix->el(i, heading_col) != 0) {
                    return i;
                }
            }
            return -1;
        }
        void solution () {
            int done_rows = 0;
            int rows = matrix -> get_rows();

            cout << "Matrix: " << endl;

            matrix -> printMatrix();

            while(done_rows < rows) {

                setHeadingElem();

                if(fabs(heading_elem) < eps) {
                    int swap_row = searhForSwapRow();

                    if(swap_row == -1) {
                        while(fabs(heading_elem) < eps) {
                            setHeadingElem();
                        }
                    } else {
                        matrix -> swap_rows(heading_row, swap_row);
                        heading_elem = matrix -> el(heading_row, heading_col);
                    }
                }

                cout << endl << "Head element is: " << heading_elem << " on " << heading_row << " row and " << "column " << heading_col << endl << endl;


                if(heading_elem != 1) {
                    float devider = heading_elem;
                    matrix -> devideRow(heading_row, devider);
                    matrix -> printMatrix();
                }

                matrix -> transformRows(heading_row, heading_col);

                removeEmptyRows();

                matrix -> printMatrix();

                cout << endl;

                done_rows++;
            }

            matrix -> printEquations();
        }

};
class Gauss: public Method{
    private:
        int* colums_done;
        int* rows_done;
        int rows_done_count;
        int cols_done_count;
        float heading_elem;
        int heading_col;
        int heading_row;
        void setHeadingElem() {
            int rows = matrix -> get_rows();
            int working_cols = matrix -> get_cols();
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < working_cols; j++) {
                    if(fabs(heading_elem) > fabs(matrix->el(i,j))
                    && !isInArray(rows_done, rows_done_count, i)
                    && !isInArray(colums_done, cols_done_count, j)
                    && matrix->el(i,j) != 0) {
                        heading_elem = matrix->el(i,j);
                        heading_col = j;
                        heading_row = i;
                    }
                }
            }
        }
        bool isInArray(int* arr, int size, int search) {
            for(int i = 0; i < size; i++) {
                if(arr[i] == search) {
                    return true;
                }
            }
            return false;
        } 
    public:
        Gauss(Matrix* m) : Method(m) {
            rows_done_count = 0;
            cols_done_count = 0;
            heading_row = 0;
            heading_col = 0;
            heading_elem = -100000;
            colums_done = new int[matrix -> get_cols()];
            rows_done = new int[matrix -> get_rows()];
        }
        ~Gauss() {
            delete [] colums_done;
            delete [] rows_done;
        }
        void solution() {

            int woking_cols = matrix -> working_cols();
            int rows = matrix -> get_rows();
            cout << "Matrix: " << endl;
            matrix -> printMatrix();

            while(cols_done_count < woking_cols && rows_done_count < rows) {

                heading_row = 0;
                heading_col = 0;
                heading_elem = -10000;

                setHeadingElem();

                cout << endl << "Head element is: " << heading_elem << " on " << heading_row << " row and " << " column " << heading_col << endl << endl;

                matrix -> transformRows(heading_row, heading_col);

                rows_done[rows_done_count] = heading_row;
                colums_done[cols_done_count] = heading_col;

                cols_done_count++;
                rows_done_count++;

                removeEmptyRows();

                matrix -> printMatrix();

                cout << endl;
            }

            matrix -> printEquations();
        }
};

int main () {
    int cols, rows;

    cout << "Enter column size: " << endl;
    cin >> cols;

    cout << "Enter row size: " << endl;
    cin >> rows;

    Matrix m = Matrix(rows, cols);

    m.readMatrix();

    int mode;
    cout << "Mode: " << endl;
    cout << " 0 - Gauss for djurkane only" << endl;
    cout << " 1 - Gauss for equations (you must enter extended matrix)" << endl;
    cout << " 2 - Gauss-Jordan for djurkane only" << endl;
    cout << " 3 - Gauss-Jordan for eqations (you must enter a extended matrix)" << endl;

    cin >> mode;

    switch(mode) {
        case 0: {
            Gauss g_o = Gauss(&m);
            m.set_expansion(false);
            g_o.solution();
            break;
        }
        case 1: {
            Gauss g = Gauss(&m);
            m.set_expansion(true);
            g.solution();
            break;
        }
        case 2: {
            GaussJordan gj_o = GaussJordan(&m);
            m.set_expansion(false);
            gj_o.solution();
            break;
        }
        case 3: {
            GaussJordan gj = GaussJordan(&m);
            m.set_expansion(true);
            gj.solution();
            break;
        }
    } 

    return 0;
}