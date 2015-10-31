/*
 * LinearSolver.cpp
 *
 *  Created on: Oct 25, 2015
 *      Author: rabbit
 */
#include "LinearSolver/LinearSolver.h"
//****************************************************************************80
LinearSolver::LinearSolver(){}
//****************************************************************************80
LinearSolver::LinearSolver(const realT atol, const realT rtol) :
    atol_(atol), rtol_(rtol){}
LinearSolver::~LinearSolver(){}
