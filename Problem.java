import ilog.concert.*;
import ilog.cplex.*;
public class Problem
{
    public IloCplex cplex;
    public IloNumVar[] probAgentOne;
    public IloNumVar[] probAgentTwo;
    public IloNumVar[] paymentsAgentOne;
    public IloNumVar[] paymentsAgentTwo;
    public Problem(IloCplex cplex, IloNumVar[] probAgentOne, IloNumVar[] probAgentTwo, IloNumVar[] paymentsAgentOne, IloNumVar[] paymentsAgentTwo)
    {
        this.cplex = cplex;
        this.probAgentOne = probAgentOne;
        this.probAgentTwo = probAgentTwo;
        this.paymentsAgentOne = paymentsAgentOne;
        this.paymentsAgentTwo = paymentsAgentTwo;
    }
}
