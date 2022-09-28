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

#pragma once

#include "matrix.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

namespace phySAT
{
    struct stp_functional_params
    {
        enum function_t
        {
            calc,
            cut,
        } strategy = calc;
    };

    class stp_impl
    {
    public:
        stp_impl(string &expression, vector<string> &tt, vector<MatrixXi> &mtxvec, stp_functional_params const &ps)
            : expression(expression), tt(tt), mtxvec(mtxvec), ps(ps)
        {
        }

        void run()
        {
            switch (ps.strategy)
            {
            case stp_functional_params::calc:
                run_calc();
                break;

            case stp_functional_params::cut:
                run_cut();
                break;
            }
        }

    private:
        void run_calc()
        {
            parser_from_expression(tt, expression);
            matrix_mapping(tt, mtxvec);
            stp_exchange_judge(tt, mtxvec);
            stp_product_judge(tt, mtxvec);
        }

        void run_cut()
        {
            parser_from_expression(tt, expression);
            matrix_mapping(tt, mtxvec);
            stp_cut(tt, mtxvec);
            stp_result(mtxvec, expression);
        }

    private:
        void parser_from_expression(vector<string> &tt, string &expression)
        {
            for (int i = 0; i < expression.size(); i++)
            {

                if ((expression[i] >= 'A' && expression[i] <= 'Z') || (expression[i] >= 'a' && expression[i] <= 'z'))
                {
                    string tmp0;
                    tmp0.push_back(expression[i]);
                    tt.push_back(tmp0);
                }
                else if (expression[i] == '!')
                {
                    string tmp1 = "MN";
                    tt.push_back(tmp1);
                }
                else if (expression[i] == '(')
                {
                    string tmp2 = "MC";
                    tt.push_back(tmp2);
                }
                else if (expression[i] == '{')
                {
                    string tmp3 = "MD";
                    tt.push_back(tmp3);
                }
                else if (expression[i] == '[')
                {
                    string tmp4 = "ME";
                    tt.push_back(tmp4);
                }
                else if (expression[i] == '/')
                {
                    string tmp5 = "MI";
                    tt.push_back(tmp5);
                }
                else if (expression[i] == '<')
                {
                    string tmp6 = "MM";
                    tt.push_back(tmp6);
                }
            }
        }

        void matrix_mapping(vector<string> &tt, vector<MatrixXi> &mtxvec)
        {
            for (int ix = 0; ix < tt.size(); ix++)
            {
                if (tt[ix] == "MN")
                {
                    MatrixXi mtxtmp(2, 2);
                    mtxtmp(0, 0) = 0;
                    mtxtmp(0, 1) = 1;
                    mtxtmp(1, 0) = 1;
                    mtxtmp(1, 1) = 0;
                    mtxvec.push_back(mtxtmp);
                }
                else if (tt[ix] == "MP")
                {
                    MatrixXi mtxtmp1(2, 2);
                    mtxtmp1(0, 0) = 1;
                    mtxtmp1(0, 1) = 0;
                    mtxtmp1(1, 0) = 0;
                    mtxtmp1(1, 1) = 1;
                    mtxvec.push_back(mtxtmp1);
                }
                else if (tt[ix] == "MC")
                {
                    MatrixXi mtxtmp2(2, 4);
                    mtxtmp2(0, 0) = 1;
                    mtxtmp2(0, 1) = 0;
                    mtxtmp2(0, 2) = 0;
                    mtxtmp2(0, 3) = 0;
                    mtxtmp2(1, 0) = 0;
                    mtxtmp2(1, 1) = 1;
                    mtxtmp2(1, 2) = 1;
                    mtxtmp2(1, 3) = 1;
                    mtxvec.push_back(mtxtmp2);
                }
                else if (tt[ix] == "MD")
                {
                    MatrixXi mtxtmp3(2, 4);
                    mtxtmp3(0, 0) = 1;
                    mtxtmp3(0, 1) = 1;
                    mtxtmp3(0, 2) = 1;
                    mtxtmp3(0, 3) = 0;
                    mtxtmp3(1, 0) = 0;
                    mtxtmp3(1, 1) = 0;
                    mtxtmp3(1, 2) = 0;
                    mtxtmp3(1, 3) = 1;
                    mtxvec.push_back(mtxtmp3);
                }
                else if (tt[ix] == "ME")
                {
                    MatrixXi mtxtmp4(2, 4);
                    mtxtmp4(0, 0) = 1;
                    mtxtmp4(0, 1) = 0;
                    mtxtmp4(0, 2) = 0;
                    mtxtmp4(0, 3) = 1;
                    mtxtmp4(1, 0) = 0;
                    mtxtmp4(1, 1) = 1;
                    mtxtmp4(1, 2) = 1;
                    mtxtmp4(1, 3) = 0;
                    mtxvec.push_back(mtxtmp4);
                }
                else if (tt[ix] == "MI")
                {
                    MatrixXi mtxtmp5(2, 4);
                    mtxtmp5(0, 0) = 1;
                    mtxtmp5(0, 1) = 0;
                    mtxtmp5(0, 2) = 1;
                    mtxtmp5(0, 3) = 1;
                    mtxtmp5(1, 0) = 0;
                    mtxtmp5(1, 1) = 1;
                    mtxtmp5(1, 2) = 0;
                    mtxtmp5(1, 3) = 0;
                    mtxvec.push_back(mtxtmp5);
                }
                else if (tt[ix] == "MM")
                {
                    MatrixXi mtxtmp6(2, 8);
                    mtxtmp6(0, 0) = 1;
                    mtxtmp6(0, 1) = 1;
                    mtxtmp6(0, 2) = 1;
                    mtxtmp6(0, 3) = 0;
                    mtxtmp6(0, 4) = 1;
                    mtxtmp6(0, 5) = 0;
                    mtxtmp6(0, 6) = 0;
                    mtxtmp6(0, 7) = 0;
                    mtxtmp6(1, 0) = 0;
                    mtxtmp6(1, 1) = 0;
                    mtxtmp6(1, 2) = 0;
                    mtxtmp6(1, 3) = 1;
                    mtxtmp6(1, 4) = 0;
                    mtxtmp6(1, 5) = 1;
                    mtxtmp6(1, 6) = 1;
                    mtxtmp6(1, 7) = 1;
                    mtxvec.push_back(mtxtmp6);
                }
                else
                {
                    MatrixXi mtxtmp7(2, 1);
                    mtxtmp7(0, 0) = 1;
                    mtxtmp7(1, 0) = 0;
                    mtxvec.push_back(mtxtmp7);
                }
            }
        }

        void stp_cut(vector<string> &tt, vector<MatrixXi> &mtxvec)
        {
            vector<MatrixXi> mtx_tmp;
            vector<string> tt_tmp;
            vector<MatrixXi> result;
            int length = tt.size();
            for (int i = 0; i < (length - 1); i++)
            {
                tt_tmp.push_back(tt[0]);
                mtx_tmp.push_back(mtxvec[0]);
                if ((tt[0] != "MN") && (tt[0] != "MC") && (tt[0] != "MD") && (tt[0] != "ME") && (tt[0] != "MI") && (tt[0] != "MR") && (tt[0] != "MW") && (tt[0] != "MM"))
                {
                    tt_tmp.pop_back();
                    mtx_tmp.pop_back();
                    stp_exchange_judge(tt_tmp, mtx_tmp);
                    stp_product_judge(tt_tmp, mtx_tmp);
                    result.push_back(mtx_tmp[0]);
                    tt_tmp.clear();
                    mtx_tmp.clear();
                }
                tt.erase(tt.begin());
                mtxvec.erase(mtxvec.begin());
            }
            mtxvec.clear();
            mtxvec.assign(result.begin(), result.end());
        }

        void stp_exchange_judge(vector<string> &tt, vector<MatrixXi> &mtxvec)
        {
            for (int i = tt.size(); i > 0; i--)
            {
                for (int j = tt.size(); j > 1; j--)
                {
                    if (((tt[j - 1] != "MN") && (tt[j - 1] != "MC") && (tt[j - 1] != "MD") && (tt[j - 1] != "ME") && (tt[j - 1] != "MI") && (tt[j - 1] != "MR") && (tt[j - 1] != "MW")) && ((tt[j - 2] != "MN") && (tt[j - 2] != "MC") && (tt[j - 2] != "MD") && (tt[j - 2] != "ME") && (tt[j - 2] != "MI") && (tt[j - 2] != "MR") && (tt[j - 2] != "MW")))
                    {
                        string tmp1_tt;
                        string tmp2_tt;
                        tmp1_tt += tt[j - 2];
                        tmp2_tt += tt[j - 1];
                        if (tmp1_tt[0] > tmp2_tt[0])
                        {
                            MatrixXi matrix_w2(4, 4);
                            matrix_w2(0, 0) = 1;
                            matrix_w2(0, 1) = 0;
                            matrix_w2(0, 2) = 0;
                            matrix_w2(0, 3) = 0;
                            matrix_w2(1, 0) = 0;
                            matrix_w2(1, 1) = 0;
                            matrix_w2(1, 2) = 1;
                            matrix_w2(1, 3) = 0;
                            matrix_w2(2, 0) = 0;
                            matrix_w2(2, 1) = 1;
                            matrix_w2(2, 2) = 0;
                            matrix_w2(2, 3) = 0;
                            matrix_w2(3, 0) = 0;
                            matrix_w2(3, 1) = 0;
                            matrix_w2(3, 2) = 0;
                            matrix_w2(3, 3) = 1;
                            string tmp3_tt = "MW";
                            string tmp4_tt = tt[j - 2];
                            tt.insert(tt.begin() + (j - 2), tt[j - 1]);
                            tt.erase(tt.begin() + (j - 1));
                            tt.insert(tt.begin() + (j - 1), tmp4_tt);
                            tt.erase(tt.begin() + j);
                            tt.insert(tt.begin() + (j - 2), tmp3_tt);
                            MatrixXi tmp_mtx = mtxvec[j - 2];
                            mtxvec.insert(mtxvec.begin() + (j - 2), mtxvec[j - 1]);
                            mtxvec.erase(mtxvec.begin() + (j - 1));
                            mtxvec.insert(mtxvec.begin() + (j - 1), tmp_mtx);
                            mtxvec.erase(mtxvec.begin() + j);
                            mtxvec.insert(mtxvec.begin() + (j - 2), matrix_w2);
                        }
                    }
                    if (((tt[j - 2] != "MN") && (tt[j - 2] != "MC") && (tt[j - 2] != "MD") && (tt[j - 2] != "ME") && (tt[j - 2] != "MI") && (tt[j - 2] != "MR") && (tt[j - 2] != "MW")) && ((tt[j - 1] == "MN") || (tt[j - 1] == "MC") || (tt[j - 1] == "MD") || (tt[j - 1] == "ME") || (tt[j - 1] == "MI") || (tt[j - 1] == "MR") || (tt[j - 1] == "MW")))
                    {
                        MatrixXi tmp1;
                        tmp1 = mtxvec[j - 2];
                        MatrixXi tmp2;
                        tmp2 = mtxvec[j - 1];
                        stpm_exchange(tmp1, tmp2);
                        string tmp_tt = tt[j - 2];
                        tt.insert(tt.begin() + (j - 2), tt[j - 1]);
                        tt.erase(tt.begin() + (j - 1));
                        tt.insert(tt.begin() + (j - 1), tmp_tt);
                        tt.erase(tt.begin() + j);
                        mtxvec.insert(mtxvec.begin() + (j - 2), tmp1);
                        mtxvec.erase(mtxvec.begin() + (j - 1));
                        mtxvec.insert(mtxvec.begin() + (j - 1), tmp2);
                        mtxvec.erase(mtxvec.begin() + j);
                    }
                    if ((tt[j - 2] == tt[j - 1]) && ((tt[j - 2] != "MN") && (tt[j - 2] != "MC") && (tt[j - 2] != "MD") && (tt[j - 2] != "ME") && (tt[j - 2] != "MI") && (tt[j - 2] != "MR") && (tt[j - 2] != "MW")))
                    {
                        MatrixXi temp_mtx(4, 2);
                        temp_mtx(0, 0) = 1;
                        temp_mtx(0, 1) = 0;
                        temp_mtx(1, 0) = 0;
                        temp_mtx(1, 1) = 0;
                        temp_mtx(2, 0) = 0;
                        temp_mtx(2, 1) = 0;
                        temp_mtx(3, 0) = 0;
                        temp_mtx(3, 1) = 1;
                        mtxvec.insert(mtxvec.begin() + (j - 2), temp_mtx);
                        mtxvec.erase(mtxvec.begin() + (j - 1));
                        tt.insert(tt.begin() + (j - 2), "MR");
                        tt.erase(tt.begin() + (j - 1));
                    }
                }
            }
        }

        void stpm_exchange(MatrixXi &matrix_f, MatrixXi &matrix_b)
        {
            MatrixXi exchange_matrix;
            MatrixXi matrix_i;
            exchange_matrix = matrix_b;
            MatrixXi matrix_tmp(2, 2);
            matrix_tmp(0, 0) = 1;
            matrix_tmp(0, 1) = 0;
            matrix_tmp(1, 0) = 0;
            matrix_tmp(1, 1) = 1;
            matrix_i = matrix_tmp;
            matrix_b = matrix_f;
            matrix_f = stp_kron_product(matrix_i, exchange_matrix);
        }

        void stp_product_judge(vector<string> &tt, vector<MatrixXi> &mtxvec)
        {
            for (auto x : tt)
            {
                cout << x << " ";
            }
            cout << endl;
            MatrixXi temp0, temp1;
            for (int ix = 1; ix < tt.size(); ix++)
            {
                if ((tt[ix] == "MW") || (tt[ix] == "MN") || (tt[ix] == "MC") || (tt[ix] == "MD") || (tt[ix] == "ME") || (tt[ix] == "MI") || (tt[ix] == "MR") || (tt[ix] == "MM"))
                {
                    temp0 = mtxvec[0];
                    temp1 = mtxvec[ix];
                    mtxvec[0] = stpm_basic_product(temp0, temp1);
                }
            }
        }

        MatrixXi stpm_basic_product(MatrixXi matrix_f, MatrixXi matrix_b)
        {
            int z;
            MatrixXi result_matrix;
            MatrixXi matrix_i1;
            int n_col = matrix_f.cols();
            int p_row = matrix_b.rows();
            if (n_col % p_row == 0)
            {
                z = n_col / p_row;
                matrix_i1 = MatrixXi::eye(z);
                MatrixXi temp;
                temp = stp_kron_product(matrix_b, matrix_i1);
                result_matrix = MatrixXi::product(matrix_f, temp);
            }
            else if (p_row % n_col == 0)
            {
                z = p_row / n_col;
                MatrixXi matrix_i2;
                matrix_i2 = MatrixXi::eye(z);
                MatrixXi temp;
                temp = stp_kron_product(matrix_f, matrix_i2);
                result_matrix = MatrixXi::product(temp, matrix_b);
            }
            return result_matrix;
        }

        MatrixXi stp_kron_product(MatrixXi matrix_f, MatrixXi matrix_b)
        {
            int m = matrix_f.rows();
            int n = matrix_f.cols();
            int p = matrix_b.rows();
            int q = matrix_b.cols();
            MatrixXi dynamic_matrix(m * p, n * q);
            for (int i = 0; i < m * p; i++)
            {
                for (int j = 0; j < n * q; j++)
                {
                    dynamic_matrix(i, j) = matrix_f(i / p, j / q) * matrix_b(i % p, j % q);
                }
            }
            return dynamic_matrix;
        }

        void stp_result(vector<MatrixXi> &mtxvec, string &expression)
        {
            string stp_result;
            cout << "STPM Result for SAT:" << endl;
            int target = 1;
            int count = 0;
            for (int i = 0; i < expression.size(); i++)
            {
                if ((expression[i] >= 'A' && expression[i] <= 'Z') || (expression[i] >= 'a' && expression[i] <= 'z'))
                {
                    count += 1;
                }
            }
            stp_result_enumeration(mtxvec, target, stp_result);
            vector<string> result;
            string tmp;
            int count_result = 0;
            for (int j = 0; j < stp_result.size(); j++)
            {
                if (stp_result[j] == '\n')
                {
                    if (tmp.size() == count)
                    {
                        result.push_back(tmp);
                        count_result += 1;
                    }
                    if (tmp.size() > count)
                    {
                        int tmp0 = tmp.size() - count - 1;
                        tmp.erase((count - tmp0), tmp0 + 1);
                        result.push_back(tmp);
                        count_result += 1;
                    }
                }
                tmp += stp_result[j];
            }
            int count_r = 0;
            for (string i : result)
            {
                cout << i << " ";
                if (count_r == 10)
                {
                    count_r = 0;
                    cout << endl;
                }
                count_r++;
            }
            cout << endl
                 << "Number of all Result:" << count_result << endl;
            if (count_result > 0)
            {
                cout << "SATISFIABLE" << endl;
            }
            else
            {
                cout << "UNSATISFIABLE" << endl;
            }
        }

        void stp_result_enumeration(vector<MatrixXi> &mtxvec, int &target, string &stp_result)
        {
            int target_tmp;
            int n = mtxvec[0].cols();
            target_tmp = target;
            for (int j = 0; j < n; j++)
            {
                if (mtxvec[0](0, j) == target_tmp)
                {
                    if (j < (n / 2))
                    {
                        string tmp = "1";
                        stp_result += tmp;
                        if (j < (n / 4))
                        {
                            target = 1;
                            if (mtxvec.size() > 1)
                            {
                                vector<MatrixXi> temp;
                                temp.assign(mtxvec.begin() + 1, mtxvec.end());
                                stp_result_enumeration(temp, target, stp_result);
                            }
                            if (mtxvec.size() == 1)
                            {
                                string tmp = "1";
                                stp_result += tmp;
                                string tmp1 = "\n";
                                stp_result += tmp1;
                            }
                        }
                        if (j >= (n / 4))
                        {
                            target = 0;
                            if (mtxvec.size() > 1)
                            {
                                vector<MatrixXi> temp;
                                temp.assign(mtxvec.begin() + 1, mtxvec.end());
                                stp_result_enumeration(temp, target, stp_result);
                            }
                            if (mtxvec.size() == 1)
                            {
                                string tmp = "0";
                                stp_result += tmp;
                                string tmp1 = "\n";
                                stp_result += tmp1;
                            }
                        }
                    }
                    if (j >= (n / 2))
                    {
                        string tmp = "0";
                        stp_result += tmp;
                        if (j < ((3 * n) / 4))
                        {
                            target = 1;
                            if (mtxvec.size() > 1)
                            {
                                vector<MatrixXi> temp;
                                temp.assign(mtxvec.begin() + 1, mtxvec.end());
                                stp_result_enumeration(temp, target, stp_result);
                            }
                            if (mtxvec.size() == 1)
                            {
                                string tmp = "1";
                                stp_result += tmp;
                                string tmp1 = "\n";
                                stp_result += tmp1;
                            }
                        }
                        if (j >= ((3 * n) / 4))
                        {
                            target = 0;
                            if (mtxvec.size() > 1)
                            {
                                vector<MatrixXi> temp;
                                temp.assign(mtxvec.begin() + 1, mtxvec.end());
                                stp_result_enumeration(temp, target, stp_result);
                            }
                            if (mtxvec.size() == 1)
                            {
                                string tmp = "0";
                                stp_result += tmp;
                                string tmp1 = "\n";
                                stp_result += tmp1;
                            }
                        }
                    }
                }
            }
        }

    private:
        string &expression;
        vector<string> &tt;
        vector<MatrixXi> &mtxvec;
        stp_functional_params const &ps;
    };

    void stp(string &expression, vector<string> &tt, vector<MatrixXi> &mtxvec, stp_functional_params const &ps = {})
    {
        stp_impl p(expression, tt, mtxvec, ps);
        p.run();
    }
}
