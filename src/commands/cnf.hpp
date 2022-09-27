/* phySAT: This is a SAT solver based on Semi-Tonser Product
 * Copyright (C) 2022 */

/**
 * @file cnf.hpp
 *
 * @brief Semi-Tensor Product based SAT solver for CNF solving
 *
 * @author Hongyang Pan
 * @since  2022/09/27
*/

#ifndef CNF_HPP
#define CNF_HPP

#include "../core/cnf.hpp"
#include <alice/alice.hpp>

using namespace std;

namespace alice
{
    class cnf_command : public command
    {
    public:
        explicit cnf_command(const environment::ptr &env) : command(env, " CNF solving ")
        {
            add_option("filename,-f", filename, "the input file name");
        }

    protected:
        void execute()
        {
            clock_t start, end;
            double totalTime;
            ifstream fin(filename);
            stringstream buffer;
            buffer << fin.rdbuf();
            string str(buffer.str());
            string expression;
            vector<int> expre;
            for (int i = 6; i < (str.size() - 1); i++)
            {
                expression += str[i];
                if ((str[i] == ' ') || (str[i] == '\n'))
                {
                    expression.pop_back();
                    int intstr = atoi(expression.c_str());
                    expre.push_back(intstr);
                    expression.clear();
                }
            }
            vector<int> &exp = expre;
            vector<string> &t = tt;
            vector<MatrixXi> &mtxvec = vec;
            int count = 0;
            start = clock();
            phySAT::stp_cnf(exp, t, mtxvec);
            end = clock();
            totalTime = (double)(end - start) / CLOCKS_PER_SEC; 
            if (t.size() > 0)
            {
                cout << "Semi-Tensor Product Result : " << endl;
                for (string i : t)
                {
                    cout << i << " ";
                    count+=1;
                    if (count == 10)
                    {
                        cout << endl;
                        count = 0;
                    }
                }
                cout << endl;
                cout << "Result number : " << t.size()<<endl; 
                cout << "SATISFIABLE" << endl;
            }
            else
            {
                cout << "UNSATISFIABLE" << endl;
            }
            cout.setf(ios::fixed);
            cout << "[CPU time]:  " << setprecision(3) << totalTime << "  seconds" << endl;
            fin.close();
        }

    private:
        string filename;
        string tmp;
        vector<int> expre;
        vector<string> tt;
        vector<MatrixXi> vec;
    };

    ALICE_ADD_COMMAND(cnf, "SAT")

}

#endif