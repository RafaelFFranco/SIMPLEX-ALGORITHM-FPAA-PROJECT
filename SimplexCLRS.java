/*
 * Algoritmo Simplex — implementação baseada no CLRS (Cap. 29)
 * Cormen, Leiserson, Rivest, Stein — "Introduction to Algorithms", 3ª ed.
 *
 * Notação do livro (Seção 29.3):
 *   N — índices das variáveis não-básicas
 *   B — índices das variáveis básicas
 *   A — coeficientes das restrições na forma de folgas
 *   b — lados direitos das restrições
 *   c — coeficientes da função objetivo
 *   v — constante acumulada do valor objetivo
 */

import java.util.ArrayList;
import java.util.List;

public class SimplexCLRS {

    enum Status {
        OPTIMAL,      // solução ótima encontrada
        INFEASIBLE,   // região viável vazia
        UNBOUNDED     // objetivo cresce indefinidamente
    }

    static class SimplexResult {
        Status status;
        double objective;
        double[] x;
    }

    static class Simplex {

        // m = número de restrições, n = número de variáveis originais
        private int m, n;

        // Forma de folgas 
        private List<Integer> N, B;
        private List<List<Double>> A;
        private List<Double> b;
        private List<Double> c;
        private double v;

        private static final double EPS = 1e-9;
        private static final double INF = 1e18;

        Simplex(double[][] Ain, double[] bin, double[] cin) {
            this.m = bin.length;
            this.n = cin.length;

            N = new ArrayList<>();
            for (int j = 0; j < n; j++) N.add(j);

            B = new ArrayList<>();
            for (int i = 0; i < m; i++) B.add(n + i);

            A = new ArrayList<>();
            for (int i = 0; i < m; i++) {
                List<Double> row = new ArrayList<>();
                for (int j = 0; j < n; j++) row.add(Ain[i][j]);
                A.add(row);
            }

            b = new ArrayList<>();
            for (double val : bin) b.add(val);

            c = new ArrayList<>();
            for (double val : cin) c.add(val);

            v = 0.0;
        }

        // PIVOT — troca a variável básica da linha l pela não-básica da coluna e
        private void pivot(int l, int e) {
            int nv = N.size();
            double coef = A.get(l).get(e);

            b.set(l, b.get(l) / coef);

            for (int j = 0; j < nv; j++)
                if (j != e)
                    A.get(l).set(j, A.get(l).get(j) / coef);

            A.get(l).set(e, 1.0 / coef);

            // propaga a troca para as outras equações
            for (int i = 0; i < m; i++) {
                if (i == l) continue;

                double aie = A.get(i).get(e);
                if (Math.abs(aie) < EPS) continue;

                b.set(i, b.get(i) - aie * b.get(l));

                for (int j = 0; j < nv; j++)
                    if (j != e)
                        A.get(i).set(j, A.get(i).get(j) - aie * A.get(l).get(j));

                A.get(i).set(e, -aie * A.get(l).get(e));
            }

            double ce = c.get(e);
            if (Math.abs(ce) > EPS) {
                v += ce * b.get(l);

                for (int j = 0; j < nv; j++)
                    if (j != e)
                        c.set(j, c.get(j) - ce * A.get(l).get(j));

                c.set(e, -ce * A.get(l).get(e));
            }

            // quem saiu da base vira não-básica, quem entrou vira básica
            int tmp = B.get(l);
            B.set(l, N.get(e));
            N.set(e, tmp);
        }

        // INITIALIZE-SIMPLEX — garante uma solução básica viável inicial,
        private boolean initializeSimplex() {

            int k = 0;
            for (int i = 1; i < m; i++)
                if (b.get(i) < b.get(k)) k = i;

            // se o menor b já é não-negativo, a origem é viável
            if (b.get(k) >= -EPS) return true;

            List<Integer> N_orig = new ArrayList<>(N);
            List<Integer> B_orig = new ArrayList<>(B);
            List<List<Double>> A_orig = deepCopy(A);
            List<Double> b_orig = new ArrayList<>(b);
            List<Double> c_orig = new ArrayList<>(c);
            double v_orig = v;

            // problema auxiliar: introduz x0 e maximiza -x0
            for (int i = 0; i < m; i++)
                A.get(i).add(1.0);
            N.add(n + m);

            int nvAux = N.size();
            c = new ArrayList<>();
            for (int j = 0; j < nvAux; j++) c.add(0.0);
            c.set(nvAux - 1, -1.0);
            v = 0.0;

            int eAuxInit = N.size() - 1;
            pivot(k, eAuxInit);

            while (true) {
                int eAux = -1;
                int bestIdx = Integer.MAX_VALUE;
                for (int j = 0; j < N.size(); j++) {
                    if (c.get(j) > EPS && N.get(j) < bestIdx) {
                        bestIdx = N.get(j);
                        eAux = j;
                    }
                }
                if (eAux == -1) break;

                int lAux = -1;
                double minRatio = INF;
                int minBIdx = Integer.MAX_VALUE;
                for (int i = 0; i < m; i++) {
                    if (A.get(i).get(eAux) > EPS) {
                        double ratio = b.get(i) / A.get(i).get(eAux);
                        if (ratio < minRatio - EPS ||
                            (Math.abs(ratio - minRatio) < EPS && B.get(i) < minBIdx)) {
                            minRatio = ratio;
                            minBIdx = B.get(i);
                            lAux = i;
                        }
                    }
                }
                if (lAux == -1) break;

                pivot(lAux, eAux);
            }

            // o problema original é inviável se o ótimo de L_aux for negativo
            if (v < -EPS) {
                N = N_orig; B = B_orig; A = A_orig;
                b = b_orig; c = c_orig; v = v_orig;
                return false;
            }

            // caso degenerado: x0 ainda básica, mas com valor 0
            for (int i = 0; i < m; i++) {
                if (B.get(i) == n + m) {
                    int eDeg = -1;
                    for (int j = 0; j < N.size(); j++) {
                        if (N.get(j) != n + m && Math.abs(A.get(i).get(j)) > EPS) {
                            eDeg = j;
                            break;
                        }
                    }
                    if (eDeg == -1) {
                        // linha redundante, sem efeito sobre a região viável
                        B.remove(i);
                        b.remove(i);
                        A.remove(i);
                        m--;
                        i--;
                    } else {
                        pivot(i, eDeg);
                    }
                    break;
                }
            }

            // remove x0 de N e de todas as colunas de A
            int colX0 = -1;
            for (int j = 0; j < N.size(); j++) {
                if (N.get(j) == n + m) { colX0 = j; break; }
            }
            if (colX0 != -1) {
                N.remove(colX0);
                c.remove(colX0);
                for (int i = 0; i < m; i++)
                    A.get(i).remove(colX0);
            }

            // reconstrói a função objetivo original sobre a base atual
            c = new ArrayList<>();
            for (int j = 0; j < N.size(); j++) c.add(0.0);
            v = v_orig;

            for (int jj = 0; jj < n; jj++) {
                double cjOrig = c_orig.get(jj);
                if (Math.abs(cjOrig) < EPS) continue;

                boolean found = false;
                for (int j = 0; j < N.size(); j++) {
                    if (N.get(j) == jj) {
                        c.set(j, c.get(j) + cjOrig);
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    // jj está em B: substitui sua equação na função objetivo
                    for (int i = 0; i < m; i++) {
                        if (B.get(i) == jj) {
                            v += cjOrig * b.get(i);
                            for (int j = 0; j < N.size(); j++)
                                c.set(j, c.get(j) - cjOrig * A.get(i).get(j));
                            break;
                        }
                    }
                }
            }

            return true;
        }

        // SIMPLEX — Fase 1 garante viabilidade, Fase 2 otimiza 
        SimplexResult solve() {
            SimplexResult result = new SimplexResult();

            if (!initializeSimplex()) {
                result.status = Status.INFEASIBLE;
                return result;
            }

            while (true) {

                // regra de Bland: entre as candidatas, escolhe a de menor índice (evita ciclagem)
                int e = -1;
                int bestIdx = Integer.MAX_VALUE;
                for (int j = 0; j < N.size(); j++) {
                    if (c.get(j) > EPS && N.get(j) < bestIdx) {
                        bestIdx = N.get(j);
                        e = j;
                    }
                }

                if (e == -1) break; // nenhuma variável melhora o objetivo: ótimo atingido

                // teste da razão mínima: define quem sai da base
                int l = -1;
                double minRatio = INF;
                int minBIdx = Integer.MAX_VALUE;

                for (int i = 0; i < m; i++) {
                    if (A.get(i).get(e) > EPS) {
                        double ratio = b.get(i) / A.get(i).get(e);
                        if (ratio < minRatio - EPS ||
                            (Math.abs(ratio - minRatio) < EPS && B.get(i) < minBIdx)) {
                            minRatio = ratio;
                            minBIdx = B.get(i);
                            l = i;
                        }
                    }
                }

                if (l == -1) {
                    result.status = Status.UNBOUNDED; // nada limita o crescimento da variável que entra
                    return result;
                }

                pivot(l, e);
            }

            result.status = Status.OPTIMAL;
            result.objective = v;
            result.x = new double[n];

            for (int i = 0; i < m; i++)
                if (B.get(i) < n)
                    result.x[B.get(i)] = b.get(i);

            return result;
        }

        private static List<List<Double>> deepCopy(List<List<Double>> src) {
            List<List<Double>> copy = new ArrayList<>();
            for (List<Double> row : src) copy.add(new ArrayList<>(row));
            return copy;
        }
    }
}