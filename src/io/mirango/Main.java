package io.mirango;

public class Main {
    // Move any terminal states to the beginning (make sure to move columns, too)
    public static int[][] moveTerminal(int[][] fuel) {
        int terminal = 0; // There are no terminal states to start (helps us order it!)
        int size = fuel.length; // Max loop

        // find terminal states and move them
        for(int i = 0; i < size; i++) {
            if(isTerminalState(fuel[i])) {
                fuel = swap(terminal, i, fuel); // if it's a terminal state swap it
                fuel[terminal][terminal] = 1; // add the identity
                terminal++;
            }
        }
        return fuel;
    }

    public static boolean isTerminalState(int[] fuel) {
        int size = fuel.length;
        for(int i = 0; i < size; i++) {
            if(fuel[i] != 0) {
                return false;
            }
        }
        return true;
    }

    // Handle swaps here
    public static int[][] swap(int terminal,int index, int[][] fuel) {
        // Row swap
        int[] swap = fuel[index];
        fuel[index] = fuel[terminal];
        fuel[terminal] = swap;

        // Column swap
        for(int i = 0; i < fuel.length; i++) {
            int temp = fuel[i][terminal];
            fuel[i][terminal] = fuel[i][index];
            fuel[i][index] = temp;
        }
        return fuel;
    }

    // Create a new int[][] for the denominators
    // Maintain the original for the numerators
    public static int[] absorbantMarkov(int[][] fuel) {
        // PREP
        // Count terminal states
        int terminal = 0;
        for(int i = 0; i < fuel.length; i++) {
            if(isTerminalState(fuel[i])) {
                terminal++;
            }
            else {
                continue;
            }
        }

        // Move terminal states and add identity to terminal states
        fuel = moveTerminal(fuel);

        // Fractional Matrix
        Fraction[][] fractionalMatrix = new Fraction[fuel.length][fuel.length];

        // Fill with the identity and zero sub-matrix
        for(int i = 0; i < terminal; i++) {
            for(int j = 0; j < fractionalMatrix.length; j++) {
                fractionalMatrix[i][j] = new Fraction(fuel[i][j]);
            }
        }
        // Fill with Fractions
        int nt = fuel.length - terminal;
        for(int i = terminal; i < fractionalMatrix.length; i++) {
            int count = 0;
            for(int j = 0; j < fuel.length; j++) {
                count += fuel[i][j];
            }
            for(int j = 0; j < fuel.length; j++) {
                fractionalMatrix[i][j] = new Fraction(fuel[i][j], count);
            }
        }

        // BEGIN ACTUAL ALGORITHM
        // I - Q
        Fraction[][] imq = new Fraction[nt][nt];
        for(int i = 0; i < nt; i++) {
            for(int j = 0; j < nt; j++) {
                imq[i][j] = Fraction.subtract(fractionalMatrix[i][j], fractionalMatrix[terminal + i][terminal + j]);
            }
        }

        // Inverse Matrix using LUP decomp
        // f = (imq)^-1
        Fraction[][] f = invert(imq);

        // Multiply Matrix
        // f * r

        // Grab r from fractionalMatrix and use in multiplication
        Fraction[][] r = new Fraction[fuel.length - terminal][terminal];
        for(int i = terminal; i < fuel.length; i++) {
            for(int j = 0; j < terminal; j++) {
                r[i - terminal][j] = Fraction.copy(fractionalMatrix[i][j]);
            }
        }

        Fraction[][] fr = multiply(f, r);

        // Return Fraction (before conversion to int array)
        Fraction[] rf = new Fraction[fr[0].length];
        // Find Common Denominator
        int s = fr[0].length;
        int lc = lcm(fr[0][0].getDenominator(), fr[0][1].getDenominator());
        for(int i = 2; i < s; i++) {
            lc = lcm(lc, fr[0][1].getDenominator());
        }

        // int array to return
        int[] result = new int[s + 1];
        for(int i = 0; i < s; i++) {
            if(fr[0][i].getNumerator() == 0) {
                continue;
            }
            result[i] = (lc / fr[0][i].getDenominator()) * fr[0][i].getNumerator();
        }
        result[s] = lc;
        return result;
    }

    public static int gcd(int a, int b) {
        return b == 0? a : gcd(b, a %b);
    }
    public static int lcm(int a, int b) {
        return ((a * b) / gcd(a,  b));
    }

    public static Fraction[][] multiply(Fraction[][]a, Fraction[][]b) {
        int ar = a.length;
        int ac = a[0].length;
        int br = b.length;
        int bc = b[0].length;

        Fraction[][] c = new Fraction[ar][bc];

        for(int i = 0; i < ar; i++) {
            for(int j = 0; j < bc; j++) {
                c[i][j] = new Fraction(0);
            }
        }

        for(int i = 0; i < ar; i++) {
            for(int j = 0; j < bc; j++) {
                for(int k = 0; k < ac; k++) {
                    c[i][j] = Fraction.add(c[i][j], Fraction.multiply(a[i][k], b[k][j]));
                }
            }
        }
        return c;
    }

    public static Fraction[][] invert(Fraction[][] a) {
        int n = a.length;
        Fraction x[][] = new Fraction[n][n];
        Fraction b[][] = new Fraction[n][n];
        int index[] = new int[n];

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++) {
                if(i == j) {
                    b[i][j] = new Fraction(1);
                }
                else {
                    b[i][j] = new Fraction(0);
                }
            }
        }

        gaussian(a, index);

        for(int i = 0; i < n - 1; ++i) {
            for(int j = i + 1; j < n; ++j) {
                for(int k = 0; k < n; ++k) {
                    b[index[j]][k] = Fraction.subtract(b[index[j]][k], Fraction.multiply(a[index[j]][i], b[index[i]][k]));
                }
            }
        }
        for(int i = 0; i < n; ++i) {
            x[n - 1][i] = Fraction.divide(b[index[n - 1]][i], a[index[n-1]][n-1]);
            for(int j = n - 2; j >= 0; --j) {
                x[j][i] = b[index[j]][i];
                for(int k = j+1; k < n; ++k) {
                    x[j][i] = Fraction.subtract(x[j][i], Fraction.multiply(a[index[j]][k], x[k][i]));
                }
                x[j][i] = Fraction.divide(x[j][i], a[index[j]][j]);
            }
        }
        return x;
    }

    public static void gaussian(Fraction[][] a, int[] index) {
        int n = index.length;
        Fraction c[] = new Fraction[n];

        for(int i = 0; i < n; i++) {
            index[i] = i;
        }

        for(int i = 0; i < n; i++) {
            Fraction c1 = new Fraction(0);
            for(int j = 0; j < n; j++) {
                Fraction c0 = Fraction.abs(a[i][j]);
                if(c0.compareTo(c1) == 1) {
                    c1 = Fraction.copy(c0);
                }
            }
            c[i] = Fraction.copy(c1);
        }

        int k = 0;
        for(int j = 0; j < n - 1; j++) {
            Fraction d = new Fraction(0);
            for(int i = j; i < n; i++) {
                Fraction d0 = Fraction.abs(a[index[i]][j]);
                d0 = Fraction.divide(d0, c[index[i]]);
                if(d0.compareTo(d) == 1) {
                    d = Fraction.copy(d0);
                    k = i;
                }
            }
            int swap = index[j];
            index[j] =  index[k];
            index[k] = swap;

            for(int i = j + 1; i < n; i++) {
                Fraction pj = Fraction.divide(a[index[i]][j], a[index[j]][j]);
                a[index[i]][j] = Fraction.copy(pj);

                for(int l = j+1; l < n; l++) {
                    a[index[i]][l] = Fraction.subtract(a[index[i]][l], Fraction.multiply(pj, a[index[j]][l]));
                }
            }
        }
    }

    public static void main(String[] args) {
        // Let's test some stuff
        //int[][] fuel = {{0,1,0,0,0,1}, {4,0,0,3,2,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
        //absorbantMarkov(fuel);
        //System.out.println(absorbantMarkov(fuel));
        int[][] f1 = {{0,2,1,0,0},{0,0,0,3,4},{0,0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
        System.out.println(absorbantMarkov(f1));
    }
}

class Fraction implements Comparable<Fraction> {
    private int numerator = 0;
    private int denominator = 1;

    public Fraction(int n) {
        numerator = n;
        denominator = 1;
    }

    public Fraction(int n, int d) {
        numerator = n;
        denominator = d;
    }

    public int getNumerator() {
        return numerator;
    }
    public int getDenominator() {
        return denominator;
    }

    public static Fraction copy(Fraction x) {
        return new Fraction(x.getNumerator(), x.getDenominator());
    }

    public static Fraction add(Fraction x, Fraction y) {
        int numerator = (x.getNumerator() * y.getDenominator()) + (x.getDenominator() * y.getNumerator());
        int denominator = x.getDenominator() * y.getDenominator();
        Fraction a = new Fraction(numerator, denominator);
        a.reduce();
        return a;
    }
    public static Fraction subtract(Fraction x, Fraction y) {
        int numerator = (x.getNumerator() * y.getDenominator()) - (x.getDenominator() * y.getNumerator());
        int denominator = x.getDenominator() * y.getDenominator();
        Fraction a = new Fraction(numerator, denominator);
        a.reduce();
        return a;
    }
    public static Fraction multiply(Fraction x, Fraction y) {
        int numerator = x.getNumerator() * y.getNumerator();
        int denominator = x.getDenominator() * y.getDenominator();
        Fraction a = new Fraction(numerator, denominator);
        a.reduce();
        return a;
    }
    // Divide is really just the reciprocal of the second number * number
    public static Fraction divide(Fraction x, Fraction y) {
        int numerator = x.getNumerator() * y.getDenominator();
        int denominator = x.getDenominator() * y.getNumerator();
        Fraction a = new Fraction(numerator, denominator);
        a.reduce();
        return a;
    }

    public static Fraction abs(Fraction x) {
        int numerator = Math.abs(x.getNumerator());
        Fraction a = new Fraction(numerator, x.getDenominator());
        return a;
    }

    // need to reduce or could overflow 32 bit int
    public void reduce() {
        int g = gcd(Math.abs(numerator), Math.abs(denominator));
        if(g == 0) {
            return;
        }
        numerator = numerator / g;
        denominator = denominator / g;
    }

    public int compareTo(Fraction x) {
        Fraction s = Fraction.subtract(this, x);
        if(numerator < 0 && x.getNumerator() > 0) {
            return -1;
        }
        if(numerator > 0 && x.getNumerator() < 0) {
            return 1;
        }
        Fraction a = Fraction.abs(this);
        Fraction b = Fraction.abs(x);
        Fraction c = Fraction.subtract(this, x);
        if(c.getNumerator() < 0) {
            return -1;
        }
        else if(c.getNumerator() > 0) {
            return 1;
        }
        else {
            return 0;
        }
    }

    public static int gcd(int a, int b) {
        return b == 0? a : gcd(b, a %b);
    }
    public static int lcm(int a, int b) {
        return ((a * b) / gcd(a,  b));
    }
}
