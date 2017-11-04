import ilog.concert.*;
import ilog.cplex.*;
public class AMDInstance {
   public static void main(String[] args) {
	//the try{}catch() is required error handling, feel free to ignore	
	try
	{
		IloCplex cplex = new IloCplex(); //Create infrastructure, you need this line
		//create vars
		IloNumVar[] vars = new IloNumVar[2]; //create an array with 2 variables

		//create two variables, '0' is lower bound, DoubleMAX_VALUE is essentially infinite 			//upper bound, Type.Float is essentially saying the variable can have a decimal 		//value, and the x1 is just giving the variable a name if you need it
		vars[0] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "x1");
		vars[1] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "x2");
		
		//create objective
		IloLinearNumExpr expr = cplex.linearNumExpr(); //NEED THIS LINE, creates an 			//expression
		expr.addTerm(3, vars[0]); // adding the term 3x1
		expr.addTerm(2, vars[1]); //adding the term 2x2
		//expression is now 3x1 + 2x2
		IloObjective objective = cplex.addMaximize(expr, "revenue");
		//adding objective function maximize 3x1 + 2x2

		//add constraints
		IloLinearNumExpr blue = cplex.linearNumExpr(); //need this for each constraint
		blue.addTerm(4, vars[0]);
		blue.addTerm(2, vars[1]);
		//4x1 + 2x2
		cplex.addLe(blue, 16, "blue"); //Le, means less than or equal to so:
		//4x1 + 2x2 <= 16, and we give this constraint the name 'blue'
		

		//we do the same things for green and red paint constraints
		IloLinearNumExpr green = cplex.linearNumExpr();
		green.addTerm(1, vars[0]);
		green.addTerm(2, vars[1]);
		cplex.addLe(green, 8, "green");
		
		IloLinearNumExpr red = cplex.linearNumExpr();
		red.addTerm(1, vars[0]);
		red.addTerm(1, vars[1]);
		cplex.addLe(red, 5, "red");

		//solve
		boolean b = cplex.solve(); //b has a value of true if solution was found
		System.out.println("Solution is " + b);

		double objectiveValue = cplex.getObjValue(); //value of objective function after 			//solution
		double variableVals [] = cplex.getValues(vars); //variableVals is an array storing 			//the values of our variables after the solution

		//Information about Solution
		IloCplex.Status status = cplex.getStatus(); //status can be "feasible, optimal, etc.
		//gives the status of solution
		System.out.println("The solution is " + status.toString());

		//Information about Objective and Variables
		//print out the objective and all the variables
		System.out.println("The objective value is " + objectiveValue);
		for (int i = 0; i < variableVals.length; i++)
		{
			System.out.println (vars[i].getName() + " = "+ variableVals[i]);
		}
		
	}
	catch(IloException e)
	{
		e.printStackTrace();
	}

   }
}
