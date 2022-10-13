/* phySAT: Semi-Tonser Product based SAT and AllSAT solver, where it can solve CNF and circuit input.
 * Copyright (C) 2022 */

/**
 * @file matrix_calc.hpp
 *
 * @brief Semi-Tensor Product of matrices calculation
 *
 * @author Hongyang Pan
 * @since  2022/09/28
 */

#ifndef CALC_HPP
#define CALC_HPP

#include "../core/matrix_calc.hpp"
#include <alice/alice.hpp>

using namespace std;

namespace alice
{
    class stp_command : public command
    {
    public:
        explicit stp_command(const environment::ptr &env) : command(env, " Semi-Tensor Product Calculation Result ")
        {
            add_option("strategy, -s", strategy, "calc = 0, cut = 1");
            add_option("filename, -f", filename, "the input txt file name");
            add_flag("--cinput,-c", "Manual input test");
            add_flag("-v,--verbose", "show statistics");
        }

    protected:
        void execute()
        {
            if (strategy == 0)
            {
                ps_stp.strategy = phySAT::stp_functional_params::calc;
            }
            else if (strategy == 1)
            {
                ps_stp.strategy = phySAT::stp_functional_params::cut;
            }
            else
            {
                assert(false);
            }

            if (is_set("cinput"))
            {
                clock_t start, end;
                double totalTime;
                cin >> expression;
                vector<string> &t = tt;
                vector<MatrixXi> &mtxvec = vec;
                string &exp = expression;
                start = clock();
                phySAT::stp(exp, t, mtxvec, ps_stp);
                end = clock();
                totalTime = (double)(end - start) / CLOCKS_PER_SEC;
                cout << "STP Result: " << endl;
                int count = 0;
                int n = mtxvec[0].cols();
                for (int p = 0; p < n; p++)
                {
                    if (mtxvec[0](0, p) == 1)
                    {
                        int i = p + 1;
                        int a[1000];
                        int b = 0;
                        int tmp = n - i;
                        while (tmp)
                        {
                            a[b] = tmp % 2;
                            tmp /= 2;
                            b++;
                        }
                        for (int q = b - 1; q >= 0; q--)
                            cout << a[q];
                        cout << " ";
                        count += 1;
                        if (count == 10)
                        {
                            cout << endl;
                            count = 0;
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
                cout << endl;
                cout.setf(ios::fixed);
                cout << "[CPU time]:  " << setprecision(3) << totalTime << "  seconds" << endl;
            }
            else
            {
                clock_t start, end;
                double totalTime;
                ifstream infile(filename);
                if (infile.is_open())
                {
                    while (!infile.eof())
                    {
                        infile >> tmp;
                        expression += tmp;
                    }
                    infile.close();
                    vector<string> &t = tt;
                    vector<MatrixXi> &mtxvec = vec;
                    string &exp = expression;
                    start = clock();
                    phySAT::stp(exp, t, mtxvec, ps_stp);
                    end = clock();
                    totalTime = (double)(end - start) / CLOCKS_PER_SEC;
                    cout.setf(ios::fixed);
                    cout << "[CPU time]:  " << setprecision(3) << totalTime << "  seconds" << endl;
                }
                else
                {
                    cerr << "Cannot open input file" << endl;
                }
            }
        }

    private:
        string filename;
        string tmp;
        string expression;
        vector<string> tt;
        vector<MatrixXi> vec;
        phySAT::stp_functional_params ps_stp;
        int strategy = 0;
    };

    ALICE_ADD_COMMAND(stp, "Matrix Calculation")

}

#endif
