// Matrix 
//  - Constructors, utilities
//  - Basic algebra (add, subtract, multiply, transpose)
//  - Scalar operations (scale, norm)
//  - Determinant via LU (partial pivoting)
//  - Row Echelon Form (REF) and Row Reduced Echelon Form (RREF)
//  - Cofactor matrix, Adjugate, Inverses (Gauss-Jordan & Adjugate)

import java.util.*;
import java.io.*;
public class Matrix {

    // initialize data fields 
    private final int rows;
    private final int cols;
    final double[][] data;  // final so they never change 
    public static final double EPS = 1e-9;

    // constructors 
    public Matrix(int rows, int cols) {
        this.rows = rows; // initilizes the arbitrary rows / columns 
        this.cols = cols;
        this.data = new double[rows][cols]; // data is the 2D array holding entries , private for direct access 
    }

    public Matrix(double[][] input) {
        this.rows = input.length;
        this.cols = input[0].length;
        this.data = new double[rows][cols]; // constructs the array , input = dimensions 
        for (int i = 0; i < rows; i++) {
            if (input[i].length != cols) {
                throw new IllegalArgumentException("All rows must have the same number of columns.");
            }
            System.arraycopy(input[i], 0, this.data[i], 0, cols); // copies each row into data 
            // makes sure that the number of rows matches columns 
        }
    }

    // Keyboard/file inputs 
    public static Matrix fromKeyboard(Scanner sc) {
        System.out.print("rows cols> ");
        int r = sc.nextInt(), c = sc.nextInt(); // reads input from scanner 
        Matrix M = new Matrix(r, c);
        System.out.println("Enter " + r + " rows of " + c + " numbers:"); // input prompt 
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                M.data[i][j] = sc.nextDouble();
        return M; // turns value into a double 
    }
    public static Matrix fromFile(String path) throws IOException { // opens file, reads and splits lines 
        java.util.List<double[]> rows = new java.util.ArrayList<>();
        int cols = -1;
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] parts = line.split("\s+");
                if (cols == -1) cols = parts.length;
                if (parts.length != cols)
                    throw new IllegalArgumentException("Inconsistent column count in file");
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++) row[j] = Double.parseDouble(parts[j]);
                rows.add(row);
            }
        }
        return new Matrix(rows.toArray(new double[0][]));
    }

    
    // utilities 
    public int getRows() { return rows; } // accessors , ease of access 
    public int getCols() { return cols; } 
    public double get(int r, int c) { return data[r][c]; }
    public void set(int r, int c, double value) { data[r][c] = value; }

    public Matrix copy() { // copies matrix contents 
        double[][] copy = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(data[i], 0, copy[i], 0, cols);
        }
        return new Matrix(copy);
    }

    // builds identity matrix 
    public static Matrix identity(int n){ Matrix I = new Matrix(n,n); 
        for(int i=0;i<n;i++) I.data[i][i]=1; return I; }

    @Override
    public String toString() { 
        StringBuilder sb = new StringBuilder(); // starts a string builder 
        for (int i = 0; i < rows; i++) {
            sb.append("[");
            for (int j = 0; j < cols; j++) {
                sb.append(String.format(Locale.US, "%8.3f", data[i][j]));
                if (j < cols - 1) sb.append(", ");
            }
            sb.append("]\n");
        }
        return sb.toString(); // prints each row with three decimal places 
    }

    // matrix algebra 
    public Matrix add(Matrix B) {
        if (rows != B.rows || cols != B.cols) // makes sure dimensions match 
            throw new IllegalArgumentException("Matrix dimensions must match for addition.");
        Matrix R = new Matrix(rows, cols); // makes new matrix 
        for (int i = 0; i < rows; i++) // iterates through each row  
            for (int j = 0; j < cols; j++) // loops through columns 
                R.data[i][j] = this.data[i][j] + B.data[i][j]; // core operation , takes element at row i , column j and adds 
        return R;
    }

    public Matrix subtract(Matrix B) {
        if (rows != B.rows || cols != B.cols)
            throw new IllegalArgumentException("Matrix dimensions must match for subtraction.");
        Matrix R = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                R.data[i][j] = this.data[i][j] - B.data[i][j]; // same thing as add 
        return R;
    }

    public Matrix multiply(Matrix B) {
        if (this.cols != B.rows) // "this" is first matrix, A 
        // B is second matrix, R is A x B 
            throw new IllegalArgumentException("Inner dimensions must match for multiplication.");
        Matrix R = new Matrix(this.rows, B.cols); // makes the result matrix take the form of rows of A x col. of B 
        // The result matrix size is determined by the outer dimensions of the two matrices being multiplied
        for (int i = 0; i < this.rows; i++) {  // loops through each row i of matrix A / controls which row of A we’re currently working on
                for (int k = 0; k < this.cols; k++) {  // loops over each column k of A / row k of B 
                double aik = this.data[i][k]; // Grabs the element from A’s i-th row and k-th column
                if (Math.abs(aik) < EPS) continue; 
                for (int j = 0; j < B.cols; j++) { // iterates across row j of B's k'th row and adds 
                    R.data[i][j] += aik * B.data[k][j]; // a[i][k] x B[k,j] into R (new matrix) [i,j]
                }
             }
        }
        return R;
    }

    public Matrix transpose() { // flips dimensions rows x cols to cols x rows 
        Matrix T = new Matrix(cols, rows);
        for (int i = 0; i < rows; i++) // loops through each row index of original matrix i = 0 first row, i = 1 second etc. 
            for (int j = 0; j < cols; j++) // same thing but with columns now 
                T.data[j][i] = data[i][j]; // swaps indices 
        return T;
    }

   
   
    

    // scalar operations
    public Matrix scale(double s) {
        Matrix R = new Matrix(rows, cols); //Creates an empty result matrix R with the same dimensions as the current matrix (this) 
        for (int i = 0; i < rows; i++) // Starts a loop over every row index i from 0 to rows - 1.
            for (int j = 0; j < cols; j++) // Inner loop over every column index j from 0 to cols - 1.
                R.data[i][j] = this.data[i][j] * s; // Multiplies the element at (i, j) by the scalar s and writes the product into the same position in R.
        return R; //Returns the new scaled matrix the original matrix is not modified 
    }

    public double norm() {
        double sum = 0.0;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                sum += data[i][j] * data[i][j]; // Squares the element at (i, j) and adds it to sum  (accumulate ∑ aᵢⱼ²).
        return Math.sqrt(sum); // Returns the square root of sum / this is the Frobenius norm ‖A‖_F = √(∑ aᵢⱼ²)
    }

   



    
    //  DETERMINANT via LU (Pivoting)  
    private static class LU { final double[][] lu; final int pivSign; LU(double[][] A, int s){ this.lu=A; this.pivSign=s; } }
    // final int pivSign stores the permutation sign (+1 or −1) from row swaps
    private static LU luDecompose(Matrix A){
        if (A.rows != A.cols) throw new IllegalArgumentException("determinant requires square matrix");
        // only square matrices 
        int n = A.rows; double[][] lu = A.toArray(); int sign = 1;
        for (int k = 0; k < n; k++) { 
            int p = k; double max = Math.abs(lu[k][k]); //best pivot row is the current row k with absolute pivot |lu[k][k]|.
            for (int i = k+1; i < n; i++){
                double v = Math.abs(lu[i][k]); if (v > max){ max = v; p = i; }
            } // //Scan rows below k to find the largest absolute entry in column k.
            // If a bigger one is found, remember its row index in p and update max. (Partial pivoting.) 
            if (max < EPS) throw new ArithmeticException("singular (zero pivot)");
            if (p != k){ double[] tmp = lu[p]; lu[p] = lu[k]; lu[k] = tmp; sign = -sign; }
            // If the best pivot is not already in row k, swap row p with row k
            //Every row swap flips the permutation sign.
            for (int i = k+1; i < n; i++){ // gauss elim 
                lu[i][k] /= lu[k][k]; double lik = lu[i][k]; // Compute the L-multiplier l_ik = lu[i][k] / lu[k][k] and store it back into lu[i][k]
                for (int j = k+1; j < n; j++) lu[i][j] -= lik * lu[k][j];//updates rest of row and forms u part 
            }
        }
        return new LU(lu, sign);
    }

    public double det(){ LU lu = luDecompose(this); double d = lu.pivSign; for (int i=0;i<rows;i++) d *= lu.lu[i][i]; return d; }
    // returns determinant 
    
    
    //Row Echelon Forms (REF & RREF)
    public Matrix ref(){
        Matrix M = copy(); int r = 0; // copy matrix M 
        for (int c = 0; c < M.cols && r < M.rows; c++){  // iterates each column and each pass tries to place 
            int piv = r; for (int i = r+1; i < M.rows; i++) // a pivot in column c at row r 
                if (Math.abs(M.data[i][c]) > Math.abs(M.data[piv][c])) piv = i;
            if (Math.abs(M.data[piv][c]) < EPS) continue; // Scan rows below r+1 to find the entry in column c with the largest absolute value
            // If a larger absolute value is found record its row index in piv  / partial pivot, and if there is no pivot, move on 
            if (piv != r){ double[] tmp = M.data[piv]; M.data[piv] = M.data[r]; M.data[r] = tmp; }
            // if the best pivot row isn’t already at r then it swaps row piv with row r to bring the pivot into position
            for (int i = r+1; i < M.rows; i++){ // eliminates values below pivot 
                double f = M.data[i][c] / M.data[r][c]; if (Math.abs(f) < EPS) continue; // if value = 0 then skip 
                M.data[i][c] = 0.0; for (int j = c+1; j < M.cols; j++) M.data[i][j] -= f * M.data[r][j];
            }
            r++;
        }
        return M;
    }

    public Matrix rref(){
        Matrix M = copy(); int r = 0; // current pivot row starts at 0 
        for (int c = 0; c < M.cols && r < M.rows; c++){ // iterate over each column and place a pivot if possible 
            int piv = -1; double max = 0.0; // initialize best pivot row to not found row , track max 
            for (int i = r; i < M.rows; i++){ double v = Math.abs(M.data[i][c]); if (v > max){ max = v; piv = i; } }
            // scan rows to find entry with max abs val 
            if (piv == -1 || max < EPS) continue;
            if (piv != r){ double[] tmp = M.data[piv]; M.data[piv] = M.data[r]; M.data[r] = tmp; }
            // if best pivot is not in place then swap it so pivot sits at (r,c) 
            double pv = M.data[r][c]; for (int j = c; j < M.cols; j++) M.data[r][j] /= pv; // read pivot value (pv)
            // divide every entry from column c to the end by pv so the pivot becomes exactly 1 
            for (int i = 0; i < M.rows; i++) if (i != r){
                double f = M.data[i][c]; if (Math.abs(f) < EPS) continue;
                //Compute the factor f we need to eliminate in column c of row i
                M.data[i][c] = 0.0; for (int j = c+1; j < M.cols; j++) M.data[i][j] -= f * M.data[r][j];
                // Zero out the column entry (i, c) explicitly 
                // Apply the row operation to the rest of the row (columns to the right) so 
                // that the pivot column becomes zero in all non-pivot rows
            }
            r++;
        }
        for (int i=0;i<M.rows;i++) for (int j=0;j<M.cols;j++) if (Math.abs(M.data[i][j]) < EPS) M.data[i][j] = 0.0;
        return M;
    }


    // Cofactors / Adjugate / Inverse (two methods) 
    public Matrix minorMatrix(int row, int col){
        if (rows != cols) throw new IllegalArgumentException("minor requires square matrix");
        Matrix M = new Matrix(rows-1, cols-1);
        for (int i=0, mi=0; i<rows; i++){
            if (i == row) continue;
            for (int j=0, mj=0; j<cols; j++){
                if (j == col) continue; M.data[mi][mj++] = this.data[i][j];
            }
            mi++;
        }
        return M;
    }

    public double cofactor(int row, int col){ double d = minorMatrix(row,col).det(); return ((row+col)%2==0) ? d : -d; }

    public Matrix cofactorMatrix(){
        if (rows != cols) throw new IllegalArgumentException("cofactor matrix requires square matrix");
        Matrix C = new Matrix(rows, cols);
        for (int i=0; i<rows; i++) for (int j=0; j<cols; j++) C.data[i][j] = cofactor(i,j);
        return C;
    }

    public Matrix adjugate(){ return cofactorMatrix().transpose(); }

    public Matrix inverseAdjugate(){
        double d = det(); if (Math.abs(d) < EPS) throw new ArithmeticException("singular matrix");
        Matrix adj = adjugate(); Matrix inv = new Matrix(rows, cols);
        for (int i=0;i<rows;i++) for (int j=0;j<cols;j++) inv.data[i][j] = adj.data[i][j] / d;
        return inv;
    }

    public Matrix inverseGaussJordan(){
        if (rows != cols) throw new IllegalArgumentException("inverse requires square matrix");
        int n = rows; Matrix aug = new Matrix(n, 2*n);
        for (int i=0;i<n;i++){ for (int j=0;j<n;j++) aug.data[i][j] = this.data[i][j]; aug.data[i][n+i] = 1.0; }
        // RREF on augmented matrix
        Matrix r = aug.rref();
        // Check left block is I
        for (int i=0;i<n;i++) for (int j=0;j<n;j++){
            double v = r.data[i][j];
            if ((i==j && Math.abs(v-1.0) > 1e-6) || (i!=j && Math.abs(v) > 1e-6))
                throw new ArithmeticException("singular matrix");
        }
        Matrix inv = new Matrix(n,n);
        for (int i=0;i<n;i++) for (int j=0;j<n;j++) inv.data[i][j] = r.data[i][n+j];
        return inv;
    }




    // Utility to dump as raw array 
    private double[][] toArray(){
        double[][] out = new double[rows][cols];
        for (int i=0;i<rows;i++) System.arraycopy(data[i], 0, out[i], 0, cols);
        return out;
    }





    public String eigenvector2x2(double lambda) {
        //only for a 2x2 matrix
        if (rows != 2 || cols != 2)
            throw new IllegalArgumentException("eigenvector2x2 requires a 2x2 matrix");

        // solve (A - λI)v = 0 for a non-zero vector v = (x, y)
        double a = data[0][0] - lambda;
        double b = data[0][1];
        double c = data[1][0];
        double d = data[1][1] - lambda;

        double x, y;

        // using the first row if it has some non-zero coefficient
        if (Math.abs(a) > Math.abs(b) && Math.abs(a) > EPS) {
            // a x + b y = 0 → x = -b/a * y, choose y = 1
            y = 1.0;
            x = -b / a;
        } else if (Math.abs(b) > EPS) {
            // a x + b y = 0 → y = -a/b * x, choose x = 1
            x = 1.0;
            y = -a / b;
        }
        // Otherwise, fall back to the second row
        else if (Math.abs(c) > Math.abs(d) && Math.abs(c) > EPS) {
            // c x + d y = 0 → x = -d/c * y, choose y = 1
            y = 1.0;
            x = -d / c;
        } else if (Math.abs(d) > EPS) {
            // c x + d y = 0 → y = -c/d * x, choose x = 1
            x = 1.0;
            y = -c / d;
        } else {
            // (A - λI) is the zero matrix - any non-zero vector is an eigenvector
            x = 1.0;
            y = 0.0;
        }

        // Scale so the largest component has magnitude 1 
        double max = Math.max(Math.abs(x), Math.abs(y));
        if (max > 0 && Math.abs(max - 1.0) > EPS) {
            x /= max;
            y /= max;
        }

        // Return a pretty string for the REPL
        return String.format(Locale.US, "(%.4f, %.4f)", x, y);
    }

    public double[] eigenvalues2x2() {
        // Only makes sense for a 2x2 matrix
        if (rows != 2 || cols != 2)
            throw new IllegalArgumentException("eigenvalues2x2 requires a 2x2 matrix");

        // A = [a b; c d]
        double a = data[0][0];
        double b = data[0][1];
        double c = data[1][0];
        double d = data[1][1];

        // Characteristic polynomial: λ^2 - (trace)λ + det = 0
        double trace = a + d;
        double det   = a * d - b * c;

        double disc = trace * trace - 4.0 * det;  // discriminant

        // If discriminant is slightly negative due to rounding, clamp to 0.
        if (disc < 0 && Math.abs(disc) <= EPS) {
            disc = 0.0;
        } else if (disc < 0) {
            // Truly complex eigenvalues (not supported in this real-only implementation)
            throw new ArithmeticException("Matrix has complex (non-real) eigenvalues");
        }

        double sqrtDisc = Math.sqrt(disc);

        // Quadratic formula: λ = (trace ± sqrt(disc)) / 2
        double lambda1 = 0.5 * (trace + sqrtDisc);
        double lambda2 = 0.5 * (trace - sqrtDisc);

        return new double[]{ lambda1, lambda2 };
    }
}
