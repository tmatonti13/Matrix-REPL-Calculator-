import java.util.*; 

class MatrixREPL {
    private static final String MENU = // menu displayed to the user with all available operations
        "\n=======================================\n" +
        "Matrix REPL\n" +
        "-----------\n" +
        "A) Make/Load matrix A\n" +
        "B) Make/Load matrix B\n" +
        "C) Print (A or B)\n" +
        "D) A + B\n" + 
        "E) A * B\n" + 
        "F) Transpose (A or B)\n" +
        "G) Determinant (A or B)\n" +
        "H) Inverse (A or B)\n" +
        "I) Row Echelon Form (A or B)\n" +
        "J) Row Canonical Form / RREF (A or B)\n" +
        "K) Cofactor & Adjugate (A or B)\n" +
        "L) Eigenvalues/Eigenvectors (2x2 exact; note for n>2)\n" +
        "Q) Quit\n" +
        "-----------\n";


    
    // two matrix variables (A and B) that the user can manipulate
    private Matrix A = null, B = null;

    // main entry point of the program
       public static void main(String[] args) { new MatrixREPL().run(); } // allows the program to run 


   // main loop, displays menu, reads user commands, and executes actions
    private void run() {
        Scanner sc = new Scanner(System.in); // scanner for reading user input
        while (true) { // infinite loop until user quits
            System.out.println(MENU);
            System.out.print("command>  ");
            String cmd = sc.next().trim().toUpperCase(Locale.US);
            try {                    
                switch (cmd) {       
                    case "A": A = makeOrLoad(sc, 'A'); break;
                    case "B": B = makeOrLoad(sc, 'B'); break;
                    case "C": printMatrix(sc); break;
                    case "D": addAB(); break;
                    case "E": multiplyAB(); break;
                    case "F": transpose(sc); break;
                    case "G": determinant(sc); break;
                    case "H": inverse(sc); break;
                    case "I": echelon(sc, false); break; // Regular Row Echelon Form
                    case "J": echelon(sc, true); break;  // Reduced Row Echelon Form
                    case "K": cofactorAdj(sc); break;
                    case "L": eigen(sc); break;
                    case "Q": System.out.println("Goodbye."); return;
                    default: System.out.println("Unknown command: " + cmd);
                } 
            } catch (Exception ex) {
                System.out.println("[error] " + ex.getMessage()); // handles errors 
            }
        }
    }


    // prompts the user to make or load a matrix from keyboard 
    private Matrix makeOrLoad(Scanner sc, char which) throws java.io.IOException {
        System.out.print("(K)eyboard or (F)ile? > ");
        String how = sc.next().trim().toUpperCase(Locale.US);
        Matrix M;
        if (how.equals("K")) M = Matrix.fromKeyboard(sc); // creates matrix from console input
        else if (how.equals("F")) { System.out.print("path> "); M = Matrix.fromFile(sc.next()); }
        else throw new IllegalArgumentException("Choose K or F");
        System.out.println("Matrix " + which + " loaded: " + M.getRows() + " x " + M.getCols());
        return M;
    }
     
    // lets the user choose which matrix (A or B) to operate on
    private Matrix pick(Scanner sc) {
        System.out.print("Pick matrix (A/B)> ");
        String t = sc.next().trim().toUpperCase(Locale.US);
        Matrix M = t.equals("A") ? A : t.equals("B") ? B : null;
        if (M == null) throw new IllegalStateException("Matrix " + t + " not loaded");
        return M;
    }

    // makes sure matrix A and B are loaded before performing operations
    private void ensureBoth() { 
        if (A == null || B == null) throw new IllegalStateException("Load both A and B first");
    }

    // just prints the matrix 
    private void printMatrix(Scanner sc) {
        System.out.println(pick(sc));
    }

    // adds both matrices and prints result 
    private void addAB() {
        ensureBoth(); 
        System.out.println("A + B =\n" + A.add(B));
    }

    // multiplies both and prints result 
    private void multiplyAB() {
        ensureBoth();
        System.out.println("A * B =\n" + A.multiply(B));
    }

    // does transpose calculation and prints output 
    private void transpose(Scanner sc) {
        System.out.println("Transpose =\n" + pick(sc).transpose());
    }

    
    private void determinant(Scanner sc) {
        Matrix M = pick(sc);
        if (M.getRows() != M.getCols()) { 
            System.out.println("[error] square required"); 
            return; 
        }
        System.out.printf(Locale.US, "det = %.6f\n", M.det());
    }

   
    private void inverse(Scanner sc) {
        Matrix M = pick(sc);
        if (M.getRows() != M.getCols()) { System.out.println("[error] square required"); return; }
        System.out.print("Method: (G)auss-Jordan or (A)djugate? > ");
        String m = sc.next().trim().toUpperCase(Locale.US);
        Matrix inv = m.equals("A") ? M.inverseAdjugate() : M.inverseGaussJordan();
        System.out.println("Inverse =\n" + inv);
    }

    private void echelon(Scanner sc, boolean reduced) {
        Matrix M = pick(sc);
        System.out.println((reduced ? "RREF" : "REF") + " =\n" + (reduced ? M.rref() : M.ref()));
    }

    private void cofactorAdj(Scanner sc) {
        Matrix M = pick(sc);
        if (M.getRows() != M.getCols()) { System.out.println("[error] square required"); return; }
        Matrix C = M.cofactorMatrix();
        System.out.println("Cofactor matrix =\n" + C);
        System.out.println("Adjugate =\n" + C.transpose());
    }

    private void eigen(Scanner sc) {
        Matrix M = pick(sc);
        if (M.getRows() != M.getCols()) { System.out.println("[error] eigen requires square matrix"); return; }
        if (M.getRows() == 2) {
            try {
                double[] ev = M.eigenvalues2x2();
                System.out.printf(Locale.US, "\u03BB1 = %.6f, \u03BB2 = %.6f\n", ev[0], ev[1]);
                System.out.println("v1 =\n" + M.eigenvector2x2(ev[0]));
                System.out.println("v2 =\n" + M.eigenvector2x2(ev[1]));
                return;
            } catch (ArithmeticException ex) {
                System.out.println("[note] 2x2 has complex/degenerate spectrum; we can add QR/power later.");
            }
        }
        System.out.println("[note] For n>2, we’ll add QR/power iteration when you’re ready.");
    }
}
