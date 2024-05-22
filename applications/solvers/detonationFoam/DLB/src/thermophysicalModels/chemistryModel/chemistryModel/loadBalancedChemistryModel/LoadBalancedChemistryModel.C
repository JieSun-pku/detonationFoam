/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/blttkgl/DLBFoam

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    LoadBalancedChemistryModel(ReactionThermo& thermo)
    : 
        StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
        balancer_(createBalancer()), 
        mapper_(createMapper(this->thermo())),
        cpuTimes_
        (
            IOobject
            (
                thermo.phasePropertyName("cellCpuTimes"),
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        ),
        refMap_
        (
            IOobject
            (
                thermo.phasePropertyName("referenceMap"),
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        )
    {
        if(balancer_.log())
        {
            cpuSolveFile_ = logFile("cpu_solve.out");
            cpuSolveFile_() << "                  time" << tab
                            << "           getProblems" << tab  
                            << "           updateState" << tab
                            << "               balance" << tab
                            << "           solveBuffer" << tab
                            << "             unbalance" << tab
                            << "               rank ID" << endl;
        }

    }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    ~LoadBalancedChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::mixtureFractionRefMapper
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createMapper
(
    const ReactionThermo& thermo
)
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                thermo.phasePropertyName("chemistryProperties"),
                thermo.db().time().constant(),
                thermo.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return mixtureFractionRefMapper(chemistryDict_tmp, thermo.composition());
}


template <class ReactionThermo, class ThermoType>
Foam::LoadBalancer
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createBalancer()
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                this->thermo().phasePropertyName("chemistryProperties"),
                this->thermo().db().time().constant(),
                this->thermo().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return LoadBalancer(chemistryDict_tmp);
}


template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    // CPU time analysis
    clockTime timer;
    scalar t_getProblems(0);
    scalar t_updateState(0);
    scalar t_balance(0);
    scalar t_solveBuffer(0);
    scalar t_unbalance(0);

    BasicChemistryModel<ReactionThermo>::correct();

    if(!this->chemistry_)
    {
        return great;
    }

    timer.timeIncrement();
    DynamicList<ChemistryProblem> allProblems = getProblems(deltaT);
    t_getProblems = timer.timeIncrement();

    RecvBuffer<ChemistrySolution> incomingSolutions;

    if(balancer_.active())
    {
        timer.timeIncrement();
        balancer_.updateState(allProblems);
        t_updateState = timer.timeIncrement();

        timer.timeIncrement();
        auto guestProblems = balancer_.balance(allProblems);
        auto ownProblems = balancer_.getRemaining(allProblems);
        t_balance = timer.timeIncrement();

        timer.timeIncrement();
        auto ownSolutions = solveList(ownProblems);
        auto guestSolutions = solveBuffer(guestProblems);
        t_solveBuffer = timer.timeIncrement();

        timer.timeIncrement();      
        incomingSolutions = balancer_.unbalance(guestSolutions);
        incomingSolutions.append(ownSolutions);
        t_unbalance = timer.timeIncrement();
    }
    else
    {
        timer.timeIncrement();
        incomingSolutions.append(solveList(allProblems));
        t_solveBuffer = timer.timeIncrement();
    }
        
    if(balancer_.log())
    {
        balancer_.printState();
        cpuSolveFile_() << setw(22)
                        << this->time().timeOutputValue()<<tab
                        << setw(22) << t_getProblems<<tab
                        << setw(22) << t_updateState<<tab
                        << setw(22) << t_balance<<tab
                        << setw(22) << t_solveBuffer<<tab
                        << setw(22) << t_unbalance<<tab
                        << setw(22) << Pstream::myProcNo()
                        << endl;
    }

    return updateReactionRates(incomingSolutions);
}


template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveSingle
(
    ChemistryProblem& problem, ChemistrySolution& solution
) const
{
    scalar timeLeft = problem.deltaT;
    scalarField c0 = problem.c;

    // Timer begins
    clockTime time;
    time.timeIncrement();


    // Calculate the chemical source terms
    while(timeLeft > small)
    {
        scalar dt = timeLeft;
        this->solve(
            problem.c,
            problem.Ti,
            problem.pi,
            dt,
            problem.deltaTChem);

        timeLeft -= dt;
    }

    solution.c_increment = (problem.c - c0) / problem.deltaT;
    solution.deltaTChem = min(problem.deltaTChem, this->deltaTChemMax_);

    // Timer ends
    solution.cpuTime = time.timeIncrement();

    solution.cellid = problem.cellid;
    solution.rhoi = problem.rhoi;
}


template <class ReactionThermo, class ThermoType>
Foam::scalar
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRates
(
    const RecvBuffer<ChemistrySolution>& solutions
)
{
    scalar deltaTMin = great;

    for(const auto& array : solutions)
    {
        for(const auto& solution : array)
        {

            for(label j = 0; j < this->nSpecie_; j++)
            {
                this->RR_[j][solution.cellid] =
                    solution.c_increment[j] * this->specieThermo_[j].W();
            }

            deltaTMin = min(solution.deltaTChem, deltaTMin);
            
            this->deltaTChem_[solution.cellid] =
                min(solution.deltaTChem, this->deltaTChemMax_);
            
            cpuTimes_[solution.cellid] = solution.cpuTime;
        }
    }

    return deltaTMin;
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveBuffer
(
    RecvBuffer<ChemistryProblem>& problems
) const
{
    // allocate the solutions buffer
    RecvBuffer<ChemistrySolution> solutions;
    
    for(auto& p : problems)
    {
        solutions.append(solveList(p));
    }
    return solutions;
}


template <class ReactionThermo, class ThermoType>
Foam::DynamicList<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveList
(
    UList<ChemistryProblem>& problems
) const
{
    DynamicList<ChemistrySolution> solutions(
        problems.size(), ChemistrySolution(this->nSpecie_));

    for(label i = 0; i < problems.size(); ++i)
    {
        solveSingle(problems[i], solutions[i]);
    }
    return solutions;
}



template <class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::DynamicList<Foam::ChemistryProblem>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getProblems
(
    const DeltaTType& deltaT
)
{
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    


    DynamicList<ChemistryProblem> solved_problems;
    DynamicList<ChemistryProblem> mapped_problems;

    solved_problems.resize(p.size(), ChemistryProblem(this->nSpecie_));

    scalarField massFraction(this->nSpecie_);
    scalarField concentration(this->nSpecie_);

    label counter = 0;
    forAll(T, celli)
    {

        if(T[celli] > this->Treact())
        {
            for(label i = 0; i < this->nSpecie_; i++)
            {
                concentration[i] = rho[celli] * this->Y_[i][celli] / this->specieThermo_[i].W();
                massFraction[i] = this->Y_[i][celli];
            }
            
            ChemistryProblem problem;
            problem.c = concentration;
            problem.Ti = T[celli];
            problem.pi = p[celli];
            problem.rhoi = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT = deltaT[celli];
            problem.cpuTime = cpuTimes_[celli];
            problem.cellid = celli;

            // This check can only be done based on the concentration as the 
            // reference temperature is not known
            if (mapper_.shouldMap(massFraction))
            {                
                mapped_problems.append(problem);
                refMap_[celli] = 1;
            }

            else 
            {
                solved_problems[counter] = problem;
                counter++;
                refMap_[celli] = 2;
            }

        }
        else
        {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }

    }

    //the real size is set here
    solved_problems.setSize(counter);

    runtime_assert(solved_problems.size() + mapped_problems.size() == p.size(), "getProblems fails");

    this->map(mapped_problems, solved_problems);
    

    return solved_problems;
}


template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::map
(
    DynamicList<ChemistryProblem>& mapped_problems, 
    DynamicList<ChemistryProblem>& solved_problems
)
{
    if (mapped_problems.size() > 0)
    {

        ChemistryProblem refProblem = mapped_problems[0];
        scalar refTemperature = refProblem.Ti;

        ChemistrySolution refSolution(this->nSpecie_);
        solveSingle(refProblem, refSolution);
        refMap_[refProblem.cellid] = 0;

        for (auto& problem : mapped_problems){

            // Check that the refmap temperature condition is also fullfilled
            if (mapper_.temperatureWithinRange(problem.Ti, refTemperature))
            {
                updateReactionRate(refSolution, problem.cellid);
                cpuTimes_[problem.cellid] = refSolution.cpuTime;
            }
            // Otherwise solve
            else 
            {
                solved_problems.append(problem);
                refMap_[problem.cellid] = 2;
            }
        }


    }
}

template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = solution.c_increment[j] * this->specieThermo_[j].W();
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}

