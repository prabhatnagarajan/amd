import ilog.concert.*;
import ilog.cplex.*;
public class Solution
{
    public IloCplex cplex;
    public boolean solve;
    public double[] probAgentOne;
    public double[] probAgentTwo;
    public double[] paymentsAgentOne;
    public double[] paymentsAgentTwo;
    public Solution(IloCplex cplex, boolean solve, double[] probAgentOne, double[] probAgentTwo, double[] paymentsAgentOne, double[] paymentsAgentTwo)
    {
        this.cplex = cplex;
        this.solve = solve;
        this.probAgentOne = probAgentOne;
        this.probAgentTwo = probAgentTwo;
        this.paymentsAgentOne = paymentsAgentOne;
        this.paymentsAgentTwo = paymentsAgentTwo;
    }
}
