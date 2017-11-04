import ilog.concert.*;
import ilog.cplex.*;
import java.util.*;
import java.io.*;
public class AutomatedMechanismDesign {

    public static class OutputData
    {
        public int numNormSamples;
        public int dirichletSamples;
        public long learnTime;
        public int numTypes1;
        public int numTypes2;
        public double correlation;
        public double confidenceLevel;
        public int violatedIRConstraints;
        public int violatedICConstraints;
        public double fullSurplus;
        public double exPostSurplus;
        public double robustSurplus;
        public double eRobustSurplus;
        public long exPostTime;
    }

    public static OutputData output = new OutputData();

   public static void main(String[] args) throws IOException{
        //System.out.println(Arrays.toString(args));
        int normalSamples = Integer.parseInt(args[2]);
        int dirichletSamples = Integer.parseInt(args[3]);
        double confidence = Double.parseDouble(args[4]);
        BayesianStatistics stats = new BayesianStatistics(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), Double.parseDouble(args[4]));
        
        output.numNormSamples = normalSamples;
        output.dirichletSamples = dirichletSamples;
        output.confidenceLevel = confidence;
        output.correlation = Double.parseDouble(args[5]);
        output.learnTime = stats.learnTime;

        double []type1 = new double[stats.agentOneTypes.length];
        double []type2 = new double[stats.agentTwoTypes.length];
        for (int i = 0; i < type1.length; i++)
        {
            type1[i] = stats.agentOneTypes[i];
        }
        for (int i = 0; i < type2.length; i++)
        {
            type2[i] = stats.agentTwoTypes[i];
        }

        output.numTypes1 = type1.length;
        output.numTypes2 = type2.length;

        double lowerBoundOne[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        double upperBoundOne[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        double lowerBoundTwo[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        double upperBoundTwo[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        double avgDist[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        double actual[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        //double upperMinusLower[][] = new double[stats.numTypesOne][stats.numTypesTwo];
        for (int i = 0; i < stats.numTypesOne; i++)
        {
          for (int k = 0; k < stats.numTypesTwo; k++)
          {
            lowerBoundOne[i][k] =  stats.sextuple.lowerOne.get(i * stats.numTypesTwo + k).getSecond();
            upperBoundOne[i][k] = stats.sextuple.upperOne.get(i * stats.numTypesTwo + k).getSecond();
            lowerBoundTwo[i][k] = stats.sextuple.lowerTwo.get(i * stats.numTypesTwo + k).getSecond();
            upperBoundTwo[i][k] = stats.sextuple.upperTwo.get(i * stats.numTypesTwo + k).getSecond();
            avgDist[i][k] = stats.sextuple.avg.get(i * stats.numTypesTwo + k).getSecond();
            //upperMinusLower[i][k] = upperBound[i][k] - lowerBound[i][k];
            actual[i][k] = stats.sextuple.actual.get(i * stats.numTypesTwo + k).getSecond();
          }
        }
        normalizeDist(avgDist);
        //normalizeDist(actual);

        /*
        TwoBidderDistributions twobid = new TwoBidderDistributions(actual, type1.length, type2.length);
        
        double epsilon = 0.025;
        
        System.out.println("PRINTING LOWER/UPPER BOUND VALUES FOR BIDDER 1 GIVEN BIDDER 2's TYPES");
        System.out.println("LOWER");
        printDoubleMatrix(getEpsilonBoundedMatrix(twobid.oneGivenTwo, epsilon, false));

        System.out.println("UPPER");
        printDoubleMatrix(getEpsilonBoundedMatrix(twobid.oneGivenTwo, epsilon, true));

        System.out.println("PRINTING LOWER/UPPER BOUND VALUES FOR BIDDER 2 GIVEN BIDDER 1's TYPES");
        System.out.println("LOWER");
        printDoubleMatrix(getEpsilonBoundedMatrix(twobid.twoGivenOne, epsilon, false));

        System.out.println("UPPER");
        printDoubleMatrix(getEpsilonBoundedMatrix(twobid.twoGivenOne, epsilon, true));

        System.out.println("PRINTING LOWER, UPPER, ACTUAL, AND ESTIMATE");
        
        System.out.println("ACTUAL");
        printDoubleMatrix(actual);

        System.out.println("AVERAGE");
        printDoubleMatrix(avgDist);

        System.out.println("LOWER");
        printDoubleMatrix(lowerBound);

        System.out.println("UPPER");
        printDoubleMatrix(upperBound);

        //System.out.println("UPPER - LOWER");
        //printDoubleMatrix(upperMinusLower);
      */
        printFullSurplus(avgDist, type1, type2);
        long before = System.currentTimeMillis();
        Solution solution = generateMechanism(type1, type2, avgDist, lowerBoundOne, upperBoundOne, lowerBoundTwo, upperBoundTwo, false);
   		  //Solution solution = generateMechanism(type1, type2, avgDist, lowerBoundOne, upperBoundOne, lowerBoundTwo, upperBoundTwo, false);
        long after = System.currentTimeMillis();
        output.exPostTime = after - before;
        //System.out.println("time("+type1.length+","+type2.length+") = " + (after - before) + " milliseconds");
        Solution solution2 = generateMechanism(type1, type2, avgDist, lowerBoundOne, upperBoundOne, lowerBoundTwo, upperBoundTwo, true);
        //Solution solution2 = generateMechanism(type1, type2, avgDist, lowerBoundOne, upperBoundOne, lowerBoundTwo, upperBoundTwo, true);
        try{
            //System.out.println("Solution 1");
            if (solution != null)
            {
                output.exPostSurplus = printSolution(type1, type2, solution);
            }
            //System.out.println();
            //System.out.println("Solution 2");
            if (solution2 != null)
            {
                output.robustSurplus = printSolution(type1, type2, solution2);
                //MUST CHANGE
                output.eRobustSurplus = output.robustSurplus * (0.95) - getWorstBNESurplus(actual, solution2);
            }
            printOutputData();
	}
        catch(IloException e)
        {
            e.printStackTrace();
        }
   }
 
   public static void printOutputData()
   {
	   System.out.print("" + output.numTypes1 + ", " + output.numTypes2 + ", " + output.numNormSamples + ", " + output.dirichletSamples + ", " + output.correlation + ", " + output.confidenceLevel + ", " + output.violatedIRConstraints + ", " + output.violatedICConstraints + ", " + output.fullSurplus + ", " + output.exPostSurplus + ", " + output.robustSurplus + ", " + output.learnTime + ", " + output.exPostTime + ", ");
   }


   public static void printFullSurplus(double dist[][], double type1[], double type2[])
   {
      //System.out.println();
      double surplus = 0;
      for (int i = 0; i < dist.length; i++)
      {
        for (int k = 0; k < dist[i].length; k++)
        {
          surplus += Math.max(type1[i], type2[k]) * dist[i][k];
        }
      }
      output.fullSurplus = surplus;
      //System.out.println("Full Surplus: " + surplus);
      //System.out.println();
   }

   //assumes:
   //single object
   //two agents
   //distribution has rows being agent one and cols being agentTwo
   public static Solution generateMechanism(double[] agentOneTypes, double[] agentTwoTypes, double[][] distribution, double [][]lowerOne, double [][]upperOne, double[][] lowerTwo, double [][] upperTwo, boolean robustIR)
   {
   		try
   		{
	   		IloCplex cplex = new IloCplex();
            cplex.setOut(null);

	   		if(!validDistribution(distribution,agentOneTypes.length, agentTwoTypes.length))
	   		{
	   			System.out.println("Bad Distribution");
	   			return null;
	   		}

          //create all the variables
	   	    //the number of variables are 2 * numagents * numtype pairs, payments and probs
  	   	  IloNumVar[] paymentsAgent1 = new IloNumVar[agentOneTypes.length * agentTwoTypes.length];
  	   	  IloNumVar[] paymentsAgent2 = new IloNumVar[agentOneTypes.length * agentTwoTypes.length];
  	   	  IloNumVar[] probAgent1 = new IloNumVar[agentOneTypes.length * agentTwoTypes.length];
          IloNumVar[] probAgent2 = new IloNumVar[agentOneTypes.length * agentTwoTypes.length];
        
          Problem problem = new Problem(cplex, probAgent1, probAgent2, paymentsAgent1, paymentsAgent2);

          double []coefficients = distributionToVector(distribution);
        
        for (int i = 0; i < agentOneTypes.length * agentTwoTypes.length; i++)
	   		{
               paymentsAgent1[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float); 
               paymentsAgent2[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
               probAgent1[i] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
               probAgent2[i] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
	   		}


            //Objective Function, maximize expcted revenue
            //create objective
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerms(coefficients, paymentsAgent1);
            expr.addTerms(coefficients, paymentsAgent2);
            IloObjective objective = cplex.addMaximize(expr, "value");

             //addExPostICAgentOne(cplex, probAgent1, paymentsAgent1, agentOneTypes, agentTwoTypes);
             //addExPostICAgentTwo(cplex, probAgent2, paymentsAgent2, agentOneTypes, agentTwoTypes);
             //addBayesianICAgentOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
             //addBayesianICAgentTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
            addSumToLeOne(cplex, probAgent1, probAgent2);
            if (robustIR)
            {
                   addExInterimIRAgentOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                   addExInterimIRAgentTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
                   addBayesianICAgentOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                   addBayesianICAgentTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
                   // addPartialInterimIROne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                   // addPartialInterimIRTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
                   // addPartialBayesianICOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                   // addPartialBayesianICTwo(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);

                   int violatedIR = addRobustIR(cplex, problem, agentOneTypes, agentTwoTypes, lowerOne, upperOne, lowerTwo, upperTwo);
                   int violatedIC = addRobustIC(cplex, problem, agentOneTypes, agentTwoTypes, lowerOne, upperOne, lowerTwo, upperTwo);
                   output.violatedIRConstraints = violatedIR;
                   output.violatedICConstraints = violatedIC;
                   //System.out.println("" + violatedIR + " Violated IR Constraints"); 
                   //System.out.println("" + violatedIC + " Violated IC Constraints"); 
            }
           else
           {
                addExPostIRAgentOne(cplex, probAgent1, paymentsAgent1, agentOneTypes, agentTwoTypes);
                addExPostIRAgentTwo(cplex, probAgent2, paymentsAgent2, agentOneTypes, agentTwoTypes);
                //addExInterimIRAgentOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                //addExInterimIRAgentTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
                //addBayesianICAgentOne(cplex, probAgent1, paymentsAgent1, distribution, agentOneTypes, agentTwoTypes);
                //addBayesianICAgentTwo(cplex, probAgent2, paymentsAgent2, distribution, agentOneTypes, agentTwoTypes);
                addExPostICAgentOne(cplex, probAgent1, paymentsAgent1, agentOneTypes, agentTwoTypes);
                addExPostICAgentTwo(cplex, probAgent2, paymentsAgent2, agentOneTypes, agentTwoTypes);
           }

           return solve(cplex, paymentsAgent1, paymentsAgent2, probAgent1, probAgent2, agentOneTypes, agentTwoTypes);

	   	}
	   	catch(IloException e)
	   	{
            e.printStackTrace();
	   	}
        return null;
   }



    public static void addPartialInterimIROne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double distribution[][],double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                double type =  agentOneTypes[i];
                double probType = probType(1, i, distribution);
                for (int j = 0; j < agentTwoTypes.length; j++)
                {
                    expectedUtility.addTerm((distribution[i][j]/probType) * type, probAgent1[i * agentTwoTypes.length + j]);
                    expectedUtility.addTerm((distribution[i][j]/probType) * (-1), paymentsAgent1[i * agentTwoTypes.length + j]);
                }
                cplex.addGe(expectedUtility, -5);
            }
    }

    public static void addPartialInterimIRTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double distribution[][],double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentTwoTypes.length; i++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                double type =  agentTwoTypes[i];
                double probType = probType(2, i, distribution);
                for (int k = 0; k < agentOneTypes.length; k++)
                {
                        expectedUtility.addTerm((distribution[k][i]/probType) * type, probAgent2[k * agentTwoTypes.length + i]);
                        expectedUtility.addTerm((distribution[k][i]/probType) * (-1), paymentsAgent2[k * agentTwoTypes.length + i]);
                }
                cplex.addGe(expectedUtility, -5);
            }

    }

    public static void addPartialExPostIROne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            //INDIVIDUAL RATIONALITY For Agent 1
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               double type = agentOneTypes[i];
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent1[i * agentOneTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent1[i * agentOneTypes.length + k]);
                  cplex.addGe(utility, -5);
               }
            }
    }

    public static void addPartialExPostIRTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
             //INDIVIDUAL RATIONALITY For Agent 2
            for (int i = 0; i < agentTwoTypes.length; i++)
            {
               double type = agentTwoTypes[i];
               for (int k = 0; k < agentOneTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent2[k * agentTwoTypes.length + i]);
                  utility.addTerm(-1, paymentsAgent2[k * agentTwoTypes.length + i]);
                  cplex.addGe(utility, -5);
               }
            }
    }

    public static void addPartialExPostICOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               double type = agentOneTypes[i];
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent1[i * agentTwoTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent1[i * agentTwoTypes.length + k]);
                  for (int j = 0; j < agentOneTypes.length; j++)
                  {
                     IloLinearNumExpr utilityLie = cplex.linearNumExpr();
                     utilityLie.addTerm(type, probAgent1[j * agentTwoTypes.length + k]);
                     utilityLie.addTerm(-1, paymentsAgent1[j * agentTwoTypes.length + k]);
                     IloLinearNumExpr constant = cplex.linearNumExpr(-5);
                     utilityLie.add(constant);     
                     cplex.addGe(utility, utilityLie);
                  }

               }
            }

    }

    public static void addPartialExPostICTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  double type = agentTwoTypes[k];
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent2[i * agentTwoTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent2[i * agentTwoTypes.length + k]);

                  for (int j = 0; j < agentTwoTypes.length; j++)
                  {
                     IloLinearNumExpr utilityLie = cplex.linearNumExpr();
                     utilityLie.addTerm(type, probAgent2[i * agentTwoTypes.length + j]);
                     utilityLie.addTerm(-1, paymentsAgent2[i * agentTwoTypes.length + j]);
                     IloLinearNumExpr constant = cplex.linearNumExpr(-5);
                     utilityLie.add(constant);      
                     cplex.addGe(utility, utilityLie);
                  }

               }
            }
    }

    public static void addPartialBayesianICOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double distribution[][], double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {       
            for (int x = 0; x < agentOneTypes.length; x++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                IloLinearNumExpr expectedUtilityLie[] = new IloLinearNumExpr[agentOneTypes.length];
                for (int y = 0; y < agentOneTypes.length; y++)
                {
                    expectedUtilityLie[y] = cplex.linearNumExpr();
                }
                double type =  agentOneTypes[x];
                double probType = probType(1, x, distribution); //check probType
                for (int i = 0; i < agentTwoTypes.length; i++)
                {
                        expectedUtility.addTerm((distribution[x][i]/probType) * type, probAgent1[x * agentTwoTypes.length + i]);
                        expectedUtility.addTerm((distribution[x][i]/probType) * (-1), paymentsAgent1[x * agentTwoTypes.length + i]);
                        for (int j = 0; j < agentOneTypes.length; j++)
                        {
                            expectedUtilityLie[j].addTerm((distribution[x][i]/probType) * type, probAgent1[j * agentTwoTypes.length + i]);
                            expectedUtilityLie[j].addTerm((distribution[x][i]/probType) * (-1), paymentsAgent1[j * agentTwoTypes.length + i]);
                        }
                }
                for (int y = 0; y < agentOneTypes.length; y++)
                {
                    IloLinearNumExpr constant = cplex.linearNumExpr(-5);
                    expectedUtilityLie[y].add(constant);                  
                    cplex.addGe(expectedUtility, expectedUtilityLie[y]);
                }
            }
    }

    public static void addPartialBayesianICTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double distribution[][], double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {       
            for (int x = 0; x < agentTwoTypes.length; x++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                IloLinearNumExpr expectedUtilityLie[] = new IloLinearNumExpr[agentTwoTypes.length];
                for (int y = 0; y < agentTwoTypes.length; y++)
                {
                    expectedUtilityLie[y] = cplex.linearNumExpr();
                }
                double type =  agentTwoTypes[x];
                double probType = probType(2, x, distribution);
                for (int i = 0; i < agentOneTypes.length; i++)
                {
                        expectedUtility.addTerm((distribution[i][x]/probType) * type, probAgent2[i * agentTwoTypes.length + x]);
                        expectedUtility.addTerm((distribution[i][x]/probType) * (-1), paymentsAgent2[i * agentTwoTypes.length + x]);                 
                        
                        for (int j = 0; j < agentTwoTypes.length; j++)
                        {
                            expectedUtilityLie[j].addTerm((distribution[i][x]/probType) * type, probAgent2[i * agentTwoTypes.length + j]);
                            expectedUtilityLie[j].addTerm((distribution[i][x]/probType) * (-1), paymentsAgent2[i * agentTwoTypes.length + j]);
                        }
                }
                
                for (int y = 0; y < agentTwoTypes.length; y++)
                {
                    IloLinearNumExpr constant = cplex.linearNumExpr(-5);
                    expectedUtilityLie[y].add(constant);       
                    cplex.addGe(expectedUtility, expectedUtilityLie[y]);
                }
            }
    }

    public static int addRobustIR(IloCplex cplex, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double[][] lowerBoundOne, double[][] upperBoundOne, double[][] lowerBoundTwo, double[][] upperBoundTwo) throws IloException
    {
        int numViolations = 0;
        boolean errorFound1, errorFound2, errorFound;
        Solution solution = null;
        do
        {
            //System.out.println("IRONE");
            solution = solve(cplex, problem.paymentsAgentOne, problem.paymentsAgentTwo, problem.probAgentOne, problem.probAgentTwo, typesAgentOne, typesAgentTwo);
            errorFound1 = identifyAndAddViolatedIRConstraintsAgentOne(solution, problem, typesAgentOne, typesAgentTwo, lowerBoundTwo, upperBoundTwo);
            //System.out.println("IRTWO");
            solution = solve(cplex, problem.paymentsAgentOne, problem.paymentsAgentTwo, problem.probAgentOne, problem.probAgentTwo, typesAgentOne, typesAgentTwo);
            errorFound2 = identifyAndAddViolatedIRConstraintsAgentTwo(solution, problem, typesAgentOne, typesAgentTwo, lowerBoundOne, upperBoundOne);
            errorFound = errorFound1 || errorFound2;
            if (errorFound)
            {
              numViolations++;
            }
        }while(errorFound);
        return numViolations;
    }

    public static int addRobustIC(IloCplex cplex, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double[][] lowerBoundOne, double[][] upperBoundOne, double[][] lowerBoundTwo, double[][] upperBoundTwo) throws IloException
    {
        int numViolations = 0;;
        boolean errorFound1, errorFound2, errorFound;
        Solution solution = null;
        do
        {
            //System.out.println("ICONE");
            solution = solve(cplex, problem.paymentsAgentOne, problem.paymentsAgentTwo, problem.probAgentOne, problem.probAgentTwo, typesAgentOne, typesAgentTwo);
            errorFound1 = identifyAndAddViolatedICConstraintsAgentOne(solution, problem, typesAgentOne, typesAgentTwo, lowerBoundTwo, upperBoundTwo);
            //System.out.println("ICTWO");
            solution = solve(cplex, problem.paymentsAgentOne, problem.paymentsAgentTwo, problem.probAgentOne, problem.probAgentTwo, typesAgentOne, typesAgentTwo);
            errorFound2 = identifyAndAddViolatedICConstraintsAgentTwo(solution, problem, typesAgentOne, typesAgentTwo, lowerBoundOne, upperBoundOne);
            errorFound = errorFound1 || errorFound2;
            if (errorFound)
            {
              numViolations++;
            }
        }while(errorFound);
        return numViolations;
    }
    
    public static boolean identifyAndAddViolatedICConstraintsAgentOne(Solution solution, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double lowerBound[][], double upperBound[][]) throws IloException
    { 
      //For all of agent one's types    
      for (int i = 0; i < typesAgentOne.length; i++)
      {
        //for each lie or alternative type of agent one
        for (int q = 0; q < typesAgentOne.length; q++)
        {
          //Create a new Cplex Object to run the linear program from paper
          IloCplex newCplex = new IloCplex();
          //Tell Cplex not to print out anything
          newCplex.setOut(null);
          //Create Variables that represent pi(theta_{-i} | theta_i)
          IloNumVar[] vars = new IloNumVar[typesAgentTwo.length];
          for (int k = 0; k < vars.length; k++)
          {
              vars[k] = newCplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
          }
          //The type of Agent One, v_i(theta_i) from paper
          double type = typesAgentOne[i];
          //create a new expression for our objective
          IloLinearNumExpr expr = newCplex.linearNumExpr();
          //create a double array of utilities, which serve as coefficients for variables in objective
          double[] utilities = new double[typesAgentTwo.length];
          for (int k = 0; k < typesAgentTwo.length; k++)
          {
              //v(theta_i) * p(theta_i, theta_{-i}) - x(theta_i, theta_{-i})
              double utility = type * solution.probAgentOne[i * typesAgentTwo.length + k] - solution.paymentsAgentOne[i * typesAgentTwo.length + k];
              //v(theta_i) * p(thetalie_i, theta_{-i}) - x(thetalie_i, theta_{-i})
              double utilityLie = type * solution.probAgentOne[q * typesAgentTwo.length + k] - solution.paymentsAgentOne[q * typesAgentTwo.length + k];
              utilities[k] = utility - utilityLie;
          }
          expr.addTerms(utilities, vars);
          IloObjective objective = newCplex.addMinimize(expr);
          for (int k = 0; k < typesAgentTwo.length; k++)
          {
              IloLinearNumExpr expression = newCplex.linearNumExpr();
              expression.addTerm(vars[k], 1);
              // double lowerSum = sumVals(1, i, lowerBound);
              // double upperSum = sumVals(1, i, upperBound);
              double lowerBoundVal = lowerBound[i][k];
              double upperBoundVal = upperBound[i][k];
              // if (upperSum == 0)
              // {
              //   lowerBoundVal = 0;
              // }
              // else
              // {
              //   lowerBoundVal = Math.max(lowerBound[i][k]/upperSum, 0);
              // }

              // if (lowerSum == 0)
              // {
              //   upperBoundVal = 1;
              // }
              // else
              // {
              //   upperBoundVal = Math.min(upperBound[i][k]/lowerSum, 1);
              // }
              newCplex.addGe(expression, lowerBoundVal);
              newCplex.addLe(expression, upperBoundVal);
          }
          IloLinearNumExpr sum = newCplex.linearNumExpr();
          double[] one = new double[vars.length];
          Arrays.fill(one, 1);
          sum.addTerms(one, vars);
          newCplex.addEq(sum, 1);
          newCplex.solve();
          double objectiveValue = newCplex.getObjValue();

          if (objectiveValue < -0.000001)
          {
              // System.out.println("violation " + objectiveValue);
              double[] varVals = newCplex.getValues(vars);
              IloLinearNumExpr utilityDiff = newCplex.linearNumExpr();
              for (int k = 0; k < typesAgentTwo.length; k++)
              {
                  utilityDiff.addTerm(varVals[k] * type, problem.probAgentOne[i * typesAgentTwo.length + k]);
                  utilityDiff.addTerm(varVals[k] * (-1), problem.paymentsAgentOne[i * typesAgentTwo.length + k]);
                  utilityDiff.addTerm(varVals[k] * type * (-1), problem.probAgentOne[q * typesAgentTwo.length + k]);
                  utilityDiff.addTerm(varVals[k], problem.paymentsAgentOne[q * typesAgentTwo.length + k]);
              }
              solution.cplex.addGe(utilityDiff, 0);
              return true;
          }
          newCplex.end();
        }
      }
      return false;
    }

    public static boolean identifyAndAddViolatedICConstraintsAgentTwo(Solution solution, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double lowerBound[][], double upperBound[][]) throws IloException
    {     
        for (int i = 0; i < typesAgentTwo.length; i++)
        {
            for (int q = 0; q < typesAgentTwo.length; q++)
            {
                IloCplex newCplex = new IloCplex();
                newCplex.setOut(null);
                IloNumVar[] vars = new IloNumVar[typesAgentOne.length];
                for (int k = 0; k < vars.length; k++)
                {
                    vars[k] = newCplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
                }
                double type = typesAgentTwo[i];
                IloLinearNumExpr expr = newCplex.linearNumExpr();
                double[] utilities = new double[typesAgentOne.length];
                for (int k = 0; k < typesAgentOne.length; k++)
                {
                    double utility = type * solution.probAgentTwo[k * typesAgentTwo.length + i] - solution.paymentsAgentTwo[k * typesAgentTwo.length + i];
                    double utilityLie = type * solution.probAgentTwo[k * typesAgentTwo.length + q] - solution.paymentsAgentTwo[k * typesAgentTwo.length + q];
                    utilities[k] = utility - utilityLie;
                }
                expr.addTerms(utilities, vars);
                IloObjective objective = newCplex.addMinimize(expr);
                for (int k = 0; k < typesAgentOne.length; k++)
                {
                    IloLinearNumExpr expression = newCplex.linearNumExpr();
                    expression.addTerm(vars[k], 1);
                    // double lowerSum = sumVals(2, i, lowerBound);
                    // double upperSum = sumVals(2, i, upperBound);
                    double lowerBoundVal = lowerBound[k][i];
                    double upperBoundVal = upperBound[k][i];
                    // if (upperSum == 0)
                    // {
                    //   lowerBoundVal = 0;
                    // }
                    // else
                    // {
                    //   lowerBoundVal = Math.max(lowerBound[k][i]/upperSum, 0);
                    // }

                    // if (lowerSum == 0)
                    // {
                    //   upperBoundVal = 1;
                    // }
                    // else
                    // {
                    //   upperBoundVal = Math.min(upperBound[k][i]/lowerSum, 1);
                    // }
                    newCplex.addGe(expression, lowerBoundVal);
                    newCplex.addLe(expression, upperBoundVal);

                    // IloLinearNumExpr lower = newCplex.linearNumExpr();
                    // IloLinearNumExpr upper = newCplex.linearNumExpr();
                    // lower.addTerm(vars[k], 1);
                    // upper.addTerm(vars[k], 1);
                    // newCplex.addGe(lower, lowerDistribution[k][i]/probTypeLower);
                    // newCplex.addLe(upper, upperDistribution[k][i]/probTypeUpper);
                }
                IloLinearNumExpr sum = newCplex.linearNumExpr();
                double[] one = new double[vars.length];
                Arrays.fill(one, 1);
                sum.addTerms(one, vars);
                newCplex.addEq(sum, 1);
                newCplex.solve();
                double objectiveValue = newCplex.getObjValue();
                if (objectiveValue < -0.000001)
                {
                    // System.out.println("violation " + objectiveValue);
                    double[] varVals = newCplex.getValues(vars);
                    IloLinearNumExpr utilityDiff = newCplex.linearNumExpr();
                    for (int k = 0; k < typesAgentOne.length; k++)
                    {
                        utilityDiff.addTerm(varVals[k] * type, problem.probAgentTwo[k * typesAgentTwo.length + i]);
                        utilityDiff.addTerm(varVals[k] * (-1), problem.paymentsAgentTwo[k * typesAgentTwo.length + i]);
                        utilityDiff.addTerm(varVals[k] * type * (-1), problem.probAgentTwo[k * typesAgentTwo.length + q]);
                        utilityDiff.addTerm(varVals[k], problem.paymentsAgentTwo[k * typesAgentTwo.length + q]);
                    }
                    solution.cplex.addGe(utilityDiff, 0);
                    return true;
                }
                newCplex.end();
            }
        }
        return false;
    }

    public static boolean identifyAndAddViolatedIRConstraintsAgentOne(Solution solution, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double lowerBound[][], double upperBound[][]) throws IloException
    {     
        for (int i = 0; i < typesAgentOne.length; i++)
        {
            IloCplex newCplex = new IloCplex();
            newCplex.setOut(null);
            IloNumVar[] vars = new IloNumVar[typesAgentTwo.length];
            for (int k = 0; k < vars.length; k++)
            {
                vars[k] = newCplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
            }
            double type = typesAgentOne[i];
            IloLinearNumExpr expr = newCplex.linearNumExpr();
            double[] utilities = new double[typesAgentTwo.length];
            for (int k = 0; k < typesAgentTwo.length; k++)
            {
                utilities[k] = type * solution.probAgentOne[i * typesAgentTwo.length + k] - solution.paymentsAgentOne[i * typesAgentTwo.length + k];
            }
            expr.addTerms(utilities, vars);
            IloObjective objective = newCplex.addMinimize(expr);
            // double probTypeLower = probType(1, i, lowerDistribution);
            // double probTypeUpper = probType(1, i, upperDistribution);
            for (int k = 0; k < typesAgentTwo.length; k++)
            {

                    IloLinearNumExpr expression = newCplex.linearNumExpr();
                    expression.addTerm(vars[k], 1);
                    // double lowerSum = sumVals(1, i, lowerBound);
                    // double upperSum = sumVals(1, i, upperBound);
                    double lowerBoundVal = lowerBound[i][k];
                    double upperBoundVal = upperBound[i][k];
                    // if (upperSum == 0)
                    // {
                    //   lowerBoundVal = 0;
                    // }
                    // else
                    // {
                    //   lowerBoundVal = Math.max(lowerBound[i][k]/upperSum, 0);
                    // }

                    // if (lowerSum == 0)
                    // {
                    //   upperBoundVal = 1;
                    // }
                    // else
                    // {
                    //   upperBoundVal = Math.min(upperBound[i][k]/lowerSum, 1);
                    // }
                    newCplex.addGe(expression, lowerBoundVal);
                    newCplex.addLe(expression, upperBoundVal);
            }
            IloLinearNumExpr sum = newCplex.linearNumExpr();
            double[] one = new double[vars.length];
            Arrays.fill(one, 1);
            sum.addTerms(one, vars);
            newCplex.addEq(sum, 1);
            newCplex.solve();
            double objectiveValue = newCplex.getObjValue();
            if (objectiveValue < -0.000001)
            {
                // System.out.println("violation " + objectiveValue);
                double[] varVals = newCplex.getValues(vars);
                IloLinearNumExpr expectedUtility = solution.cplex.linearNumExpr();
                for (int k = 0; k < typesAgentTwo.length; k++)
                {
                    expectedUtility.addTerm(varVals[k] * type, problem.probAgentOne[i * typesAgentTwo.length + k]);
                    expectedUtility.addTerm(varVals[k] * (-1), problem.paymentsAgentOne[i * typesAgentTwo.length + k]);                   
                }
                solution.cplex.addGe(expectedUtility, 0);
                return true;
            }
            newCplex.end();
        }
        return false;
    }

    public static boolean identifyAndAddViolatedIRConstraintsAgentTwo(Solution solution, Problem problem, double[] typesAgentOne, double[] typesAgentTwo, double lowerBound[][], double upperBound[][]) throws IloException
    {     
        for (int i = 0; i < typesAgentTwo.length; i++)
        {
            IloCplex newCplex = new IloCplex();
            newCplex.setOut(null);
            IloNumVar[] vars = new IloNumVar[typesAgentOne.length];
            for (int k = 0; k < vars.length; k++)
            {
                vars[k] = newCplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
            }
            double type = typesAgentTwo[i];
            IloLinearNumExpr expr = newCplex.linearNumExpr();
            double[] utilities = new double[typesAgentOne.length];
            for (int k = 0; k < typesAgentOne.length; k++)
            {
                utilities[k] = type * solution.probAgentTwo[k * typesAgentTwo.length + i] - solution.paymentsAgentTwo[k * typesAgentTwo.length+ i];
            }
            expr.addTerms(utilities, vars);
            IloObjective objective = newCplex.addMinimize(expr);
            for (int k = 0; k < typesAgentOne.length; k++)
            {
                IloLinearNumExpr expression = newCplex.linearNumExpr();
                expression.addTerm(vars[k], 1);
                // double lowerSum = sumVals(2, i, lowerBound);
                // double upperSum = sumVals(2, i, upperBound);
                double lowerBoundVal = lowerBound[k][i];
                double upperBoundVal = upperBound[k][i];
                // if (upperSum == 0)
                // {
                //   lowerBoundVal = 0;
                // }
                // else
                // {
                //   lowerBoundVal = Math.max(lowerBound[k][i]/upperSum, 0);
                // }

                // if (lowerSum == 0)
                // {
                //   upperBoundVal = 1;
                // }
                // else
                // {
                //   upperBoundVal = Math.min(upperBound[k][i]/lowerSum, 1);
                // }
                newCplex.addGe(expression, lowerBoundVal);
                newCplex.addLe(expression, upperBoundVal);
            }
            IloLinearNumExpr sum = newCplex.linearNumExpr();
            double[] one = new double[vars.length];
            Arrays.fill(one, 1);
            sum.addTerms(one, vars);
            newCplex.addEq(sum, 1);
            newCplex.solve();
            double objectiveValue = newCplex.getObjValue();
            if (objectiveValue < -0.000001)
            {
                // System.out.println("violation " + objectiveValue);
                double[] varVals = newCplex.getValues(vars);
                IloLinearNumExpr expectedUtility = solution.cplex.linearNumExpr();
                for (int k = 0; k < typesAgentOne.length; k++)
                {
                    expectedUtility.addTerm(varVals[k] * type, problem.probAgentTwo[k * typesAgentTwo.length + i]);
                    expectedUtility.addTerm(varVals[k] * (-1), problem.paymentsAgentTwo[k * typesAgentTwo.length + i]);                   
                }
                solution.cplex.addGe(expectedUtility, 0);
                return true;
            }
            newCplex.end();
        }
        return false;
    }


    public static void addSumToLeOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] probAgent2) throws IloException
    {
            //CONSTRAINTS:
            //probabilities sum to less than 1 or to 1
            for (int i = 0; i < probAgent1.length; i++)
            {
               IloLinearNumExpr probSum = cplex.linearNumExpr();
               probSum.addTerm(probAgent1[i], 1);
               probSum.addTerm(probAgent2[i], 1);
               cplex.addLe(probSum, 1);
            }

    }

    public static void addBayesianICAgentOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double distribution[][], double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {       
            for (int x = 0; x < agentOneTypes.length; x++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                IloLinearNumExpr expectedUtilityLie[] = new IloLinearNumExpr[agentOneTypes.length];
                for (int y = 0; y < agentOneTypes.length; y++)
                {
                    expectedUtilityLie[y] = cplex.linearNumExpr();
                }
                double type =  agentOneTypes[x];
                double probType = probType(1, x, distribution); //check probType
                for (int i = 0; i < agentTwoTypes.length; i++)
                {
                        expectedUtility.addTerm((distribution[x][i]/probType) * type, probAgent1[x * agentTwoTypes.length + i]);
                        expectedUtility.addTerm((distribution[x][i]/probType) * (-1), paymentsAgent1[x * agentTwoTypes.length + i]);                 
                        
                        for (int j = 0; j < agentOneTypes.length; j++)
                        {
                            expectedUtilityLie[j].addTerm((distribution[x][i]/probType) * type, probAgent1[j * agentTwoTypes.length + i]);
                            expectedUtilityLie[j].addTerm((distribution[x][i]/probType) * (-1), paymentsAgent1[j * agentTwoTypes.length + i]);
                        }
                }
                for (int y = 0; y < agentOneTypes.length; y++)
                {
                    cplex.addGe(expectedUtility, expectedUtilityLie[y]);
                }
            }
    }

    public static void addBayesianICAgentTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double distribution[][], double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {       
            for (int x = 0; x < agentTwoTypes.length; x++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                IloLinearNumExpr expectedUtilityLie[] = new IloLinearNumExpr[agentTwoTypes.length];
                for (int y = 0; y < agentTwoTypes.length; y++)
                {
                    expectedUtilityLie[y] = cplex.linearNumExpr();
                }
                double type =  agentTwoTypes[x];
                double probType = probType(2, x, distribution);
                for (int i = 0; i < agentOneTypes.length; i++)
                {
                        expectedUtility.addTerm((distribution[i][x]/probType) * type, probAgent2[i * agentTwoTypes.length + x]);
                        expectedUtility.addTerm((distribution[i][x]/probType) * (-1), paymentsAgent2[i * agentTwoTypes.length + x]);                 
                        
                        for (int j = 0; j < agentTwoTypes.length; j++)
                        {
                            expectedUtilityLie[j].addTerm((distribution[i][x]/probType) * type, probAgent2[i * agentTwoTypes.length + j]);
                            expectedUtilityLie[j].addTerm((distribution[i][x]/probType) * (-1), paymentsAgent2[i * agentTwoTypes.length + j]);
                        }
                }
                
                for (int y = 0; y < agentTwoTypes.length; y++)
                {
                    cplex.addGe(expectedUtility, expectedUtilityLie[y]);
                }
            }
    }

    public static void addExInterimIRAgentOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double distribution[][],double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                double type =  agentOneTypes[i];
                double probType = probType(1, i, distribution);
                for (int j = 0; j < agentTwoTypes.length; j++)
                {
                    expectedUtility.addTerm((distribution[i][j]/probType) * type, probAgent1[i * agentTwoTypes.length + j]);
                    expectedUtility.addTerm((distribution[i][j]/probType) * (-1), paymentsAgent1[i * agentTwoTypes.length + j]);
                }
                cplex.addGe(expectedUtility, 0);
            }
    }

    public static void addExInterimIRAgentTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double distribution[][],double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentTwoTypes.length; i++)
            {
                IloLinearNumExpr expectedUtility = cplex.linearNumExpr();
                double type =  agentTwoTypes[i];
                double probType = probType(2, i, distribution);
                for (int k = 0; k < agentOneTypes.length; k++)
                {
                        expectedUtility.addTerm((distribution[k][i]/probType) * type, probAgent2[k * agentTwoTypes.length + i]);
                        expectedUtility.addTerm((distribution[k][i]/probType) * (-1), paymentsAgent2[k * agentTwoTypes.length + i]);
                }
                cplex.addGe(expectedUtility, 0);
            }

    }


    public static double [][] getCondDist(double dist[][], double probType, int agentNum, int typeIndex)
    {
        double [][] condDist = new double [dist.length][dist[0].length];
        if (agentNum == 1)
        {
            for (int i = 0; i < dist[typeIndex].length; i++)
            {
                    condDist[typeIndex][i] = dist[typeIndex][i]/probType;
            }
            return condDist;
        }
        else if (agentNum == 2)
        {
            for (int i = 0; i < dist.length; i++)
            {
                    condDist[i][typeIndex] = dist[i][typeIndex]/probType;
            }
            return condDist; 
        }
        return condDist;
    }

    public static double sumVals(int agentNum, int index, double matrix[][])
    {
        double sum = 0;
        if (agentNum == 1)
        {
            for (int i = 0; i < matrix[index].length; i++)
            {
                sum += matrix[index][i];
            }
            return sum;
        }
        else if (agentNum == 2)
        {
            for (int i = 0; i < matrix.length; i++)
            {
                sum += matrix[i][index];
            }
            return sum;
        }
        return -1000000;
    }

    public static double probType(int agentNum, int index, double distribution[][])
    {
        double probSum = 0;
        if (agentNum == 1)
        {
            for (int i = 0; i < distribution[index].length; i++)
            {
                probSum += distribution[index][i];
            }
            return probSum;
        }
        else if (agentNum == 2)
        {
            for (int i = 0; i < distribution.length; i++)
            {
                probSum += distribution[i][index];
            }
            return probSum;
        }
        return -1000000;
    }
    public static void addExPostICAgentOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               double type = agentOneTypes[i];
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent1[i * agentTwoTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent1[i * agentTwoTypes.length + k]);
                  for (int j = 0; j < agentOneTypes.length; j++)
                  {
                     IloLinearNumExpr utilityLie = cplex.linearNumExpr();
                     utilityLie.addTerm(type, probAgent1[j * agentTwoTypes.length + k]);
                     utilityLie.addTerm(-1, paymentsAgent1[j * agentTwoTypes.length + k]);
                     cplex.addGe(utility, utilityLie);
                  }

               }
            }

    }

    public static void addExPostICAgentTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  double type = agentTwoTypes[k];
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent2[i * agentTwoTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent2[i * agentTwoTypes.length + k]);

                  for (int j = 0; j < agentTwoTypes.length; j++)
                  {
                     IloLinearNumExpr utilityLie = cplex.linearNumExpr();
                     utilityLie.addTerm(type, probAgent2[i * agentTwoTypes.length + j]);
                     utilityLie.addTerm(-1, paymentsAgent2[i * agentTwoTypes.length + j]);
                     cplex.addGe(utility, utilityLie);
                  }
               }
            }
    }


    public static void addExPostIRAgentOne(IloCplex cplex, IloNumVar[] probAgent1, IloNumVar[] paymentsAgent1, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
            //INDIVIDUAL RATIONALITY For Agent 1
            for (int i = 0; i < agentOneTypes.length; i++)
            {
               double type = agentOneTypes[i];
               for (int k = 0; k < agentTwoTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent1[i * agentTwoTypes.length + k]);
                  utility.addTerm(-1, paymentsAgent1[i * agentTwoTypes.length + k]);
                  cplex.addGe(utility, 0);
               }
            }
    }

    public static void addExPostIRAgentTwo(IloCplex cplex, IloNumVar[] probAgent2, IloNumVar[] paymentsAgent2, double agentOneTypes[], double agentTwoTypes[]) throws IloException
    {
             //INDIVIDUAL RATIONALITY For Agent 2
            for (int i = 0; i < agentTwoTypes.length; i++)
            {
               double type = agentTwoTypes[i];
               for (int k = 0; k < agentOneTypes.length; k++)
               {
                  IloLinearNumExpr utility = cplex.linearNumExpr();
                  utility.addTerm(type, probAgent2[k * agentTwoTypes.length + i]);
                  utility.addTerm(-1, paymentsAgent2[k * agentTwoTypes.length + i]);
                  cplex.addGe(utility, 0);
               }
            }
    }


   public static Solution solve(IloCplex cplex, IloNumVar[] paymentsAgent1, IloNumVar[] paymentsAgent2, IloNumVar[] probAgent1, IloNumVar[] probAgent2, double typesAgent1[], double typesAgent2[]) throws IloException
   {
            //solve
            boolean b = cplex.solve();
            //System.out.println("Solution is " + b);

            double agentOnePayments[] = cplex.getValues(paymentsAgent1);
            double agentTwoPayments[] = cplex.getValues(paymentsAgent2);
            double probabilitiesAgent1[] = cplex.getValues(probAgent1);
            double probabilitiesAgent2[] = cplex.getValues(probAgent2);

            //Information about Solution
            //IloCplex.Status status = cplex.getStatus();
            //System.out.println("The solution is " + status.toString());

            //System.out.println("The objective value is " + objectiveValue);

            //printSolution(typesAgent1, typesAgent2, agentOnePayments, agentTwoPayments, probabilitiesAgent1, probabilitiesAgent2);
            return new Solution(cplex, b, probabilitiesAgent1, probabilitiesAgent2, agentOnePayments, agentTwoPayments);   
   }

   public static double printSolution(double typesAgent1[], double typesAgent2[], Solution solution) throws IloException
   {

      //System.out.println("Solution is " + solution.solve);
      double objectiveValue = solution.cplex.getObjValue();
      //System.out.println("vars is " + (solution.cplex == null));

      //information about solution
      //System.out.println("The solution is " + solution.cplex.getStatus().toString());
      

      //System.out.print("the objective value is " + objectiveValue);


    
    /*
      System.out.println();
      System.out.println("Agent 1:");
      printPaymentsAgentOne(typesAgent1, typesAgent2, solution);
      printProbsAgentOne(typesAgent1, typesAgent2, solution);
      

      System.out.println();
      System.out.println("Agent 2:");
      System.out.println("Payments:");
      printPaymentsAgentTwo(typesAgent1, typesAgent2, solution);

      printProbsAgentTwo(typesAgent1, typesAgent2, solution);
    */
    return objectiveValue;
   }

   public static void printPaymentsAgentOne(double []typesAgent1, double []typesAgent2, Solution solution)
   {
      System.out.println("Payments:");
      for (int i = 0; i < typesAgent1.length; i++)
      {
         for (int k = 0; k < typesAgent2.length; k++)
         {
            System.out.println("paymentAgent1("+typesAgent1[i]+", "+typesAgent2[k]+") = "+ solution.paymentsAgentOne[i * typesAgent2.length + k]);
         }
      }
   }

    public static void printPaymentsAgentTwo(double []typesAgent1, double []typesAgent2, Solution solution)
    {
      System.out.println("Payments:");
      for (int i = 0; i < typesAgent1.length; i++)
      {
         for (int k = 0; k < typesAgent2.length; k++)
         {
            System.out.println("paymentAgent2("+typesAgent1[i]+", "+typesAgent2[k]+") = "+ solution.paymentsAgentTwo[i * typesAgent2.length + k]);
         }
      }
    }

    public static void printProbsAgentOne(double []typesAgent1, double []typesAgent2, Solution solution)
    {
      System.out.println("Probability that Agent 1 gets the item:");
      for (int i = 0; i < typesAgent1.length; i++)
      {
         for (int k = 0; k < typesAgent2.length; k++)
         {
            System.out.println("P("+typesAgent1[i]+", "+typesAgent2[k]+") = "+ solution.probAgentOne[i * typesAgent2.length + k]);
         }
      }
    } 

    public static void printProbsAgentTwo(double []typesAgent1, double []typesAgent2, Solution solution)
    {
      System.out.println("Probability that Agent 2 gets the item:");
      for (int i = 0; i < typesAgent1.length; i++)
      {
         for (int k = 0; k < typesAgent2.length; k++)
         {
            System.out.println("P("+typesAgent1[i]+", "+typesAgent2[k]+") = "+ solution.probAgentTwo[i * typesAgent2.length + k]);
         }
      }
    } 



   public static double[] distributionToVector(double[][] distribution)
   {
      double []vector = new double[distribution.length * distribution[0].length];
      for (int i = 0; i < distribution.length; i++)
      {
         for (int  k = 0; k < distribution[i].length; k++)
         {
            vector[i * distribution[i].length + k] = distribution[i][k ];
         }
      }
      return vector;
   }

   public static boolean validDistribution(double distribution[][], int agentOneNumTypes, int agentTwoNumTypes)
   {
   		if(distribution.length != agentOneNumTypes || distribution[0].length != agentTwoNumTypes)
   		{
   			return false;
   		}

   		for (int i = 0; i < distribution.length; i++)
   		{
   			for (int  k = 0; k < distribution[0].length; k++)
   			{
   				if(distribution[i][k] < 0)
   				{
   					return false;
   				}
   			}
   		}
        return true;
   }

   public static void normalizeDist(double distribution[][])
   {
     	double sum = 0;
   		for (int i = 0; i < distribution.length; i++)
   		{
   			for (int  k = 0; k < distribution[i].length; k++)
   			{
   				sum += distribution[i][k];
   			}
   		}

   		for (int i = 0; i < distribution.length; i++)
   		{
   			for (int  k = 0; k < distribution[i].length; k++)
   			{
   				distribution[i][k] = distribution[i][k]/sum;
   			}
   		}
   }

  public static void printDoubleMatrix(double matrix[][])
  {
    for (int i = 0; i < matrix.length; i++)
    {
      System.out.println(Arrays.toString(matrix[i]));
      System.out.println();
    }
  }

  public static double[][] getEpsilonBoundedMatrix(double dist[][], double epsilon, boolean upperBound)
  {
    double newDist[][] = new double[dist.length][dist[0].length];
    for (int i = 0; i < dist.length; i++)
    {
        for (int k = 0; k < dist[0].length; k++)
        {
            if (upperBound)
            {
                newDist[i][k] = Math.min(dist[i][k] + epsilon, 1);
            }
            else
            {
                newDist[i][k] = Math.max(dist[i][k] - epsilon, 0);
            }
        }
    }
    return newDist;
  }    
   
  //note: the distribution must be the actual distribution 
  public static double getWorstBNESurplus(double dist[][], Solution solution) throws IloException
  {
    System.out.println("I NEED SLEEP");
    IloCplex cplex = new IloCplex();
    cplex.setOut(null);
    int numTypesOne = dist.length;
    int numTypesTwo = dist[0].length;
    //create variables
    IloNumVar[] bidderOneProbs = new IloNumVar[numTypesOne * numTypesOne];
    IloNumVar[] bidderTwoProbs = new IloNumVar[numTypesTwo * numTypesTwo];
    IloNumVar[] bidderOneDecisions = new IloNumVar[numTypesOne * numTypesOne];
    IloNumVar[] bidderTwoDecisions = new IloNumVar[numTypesTwo * numTypesTwo];
    IloNumVar[] bidderOneEUs = new IloNumVar[numTypesOne * numTypesOne];
    IloNumVar[] bidderTwoEUs = new IloNumVar[numTypesTwo * numTypesTwo];
    IloNumVar[] typeUtilsOne = new IloNumVar[numTypesOne]; 
    IloNumVar[] typeUtilsTwo = new IloNumVar[numTypesTwo];
    
    for (int i = 0; i < numTypesOne * numTypesOne; i++)
    {
        bidderOneProbs[i] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
        bidderOneDecisions[i] = cplex.intVar(0, 1);
        bidderOneEUs[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
    }
    for (int i = 0; i < numTypesTwo * numTypesTwo; i++)
    {
        bidderTwoProbs[i] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
        bidderTwoDecisions[i] = cplex.intVar(0, 1);
        bidderTwoEUs[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
    }
    for (int i = 0; i < typeUtilsOne.length; i++)
    {
        typeUtilsOne[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
    }
    for (int i = 0; i < typeUtilsTwo.length; i++)
    {
        typeUtilsTwo[i] = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE, IloNumVarType.Float);
    }
    //expression (26) and (27) from paper
    for (int i = 0; i < numTypesOne; i++)
    {
        IloLinearNumExpr decOneSum = cplex.linearNumExpr();
        IloLinearNumExpr probOneSum = cplex.linearNumExpr();
        for (int k = 0; k < numTypesOne; k++)
        {
            decOneSum.addTerm(1, bidderOneDecisions[i * numTypesOne + k]);
            probOneSum.addTerm(1, bidderOneProbs[i * numTypesOne + k]);
        }
        cplex.addGe(decOneSum, 1);
        cplex.addEq(probOneSum, 1);
    }    
    for (int i = 0; i < numTypesTwo; i++)
    {
        IloLinearNumExpr decTwoSum = cplex.linearNumExpr();
        IloLinearNumExpr probTwoSum = cplex.linearNumExpr();
        for (int k = 0; k < numTypesTwo; k++)
        {
            decTwoSum.addTerm(1, bidderTwoDecisions[i * numTypesTwo + k]);
            probTwoSum.addTerm(1, bidderTwoProbs[i * numTypesTwo + k]);
        }
        cplex.addGe(decTwoSum, 1);
        cplex.addEq(probTwoSum, 1);
    }
    //expression (25) 
    for (int i = 0; i < numTypesOne * numTypesOne; i++)
    {
        IloLinearNumExpr probTerm = cplex.linearNumExpr();
        probTerm.addTerm(1, bidderOneProbs[i]);
        IloLinearNumExpr decTerm = cplex.linearNumExpr();
        decTerm.addTerm(1, bidderOneDecisions[i]);
        cplex.addLe(probTerm, decTerm);
    }
    for (int i = 0; i < numTypesTwo * numTypesTwo; i++)
    {
        IloLinearNumExpr probTerm = cplex.linearNumExpr();
        probTerm.addTerm(1, bidderTwoProbs[i]);
        IloLinearNumExpr decTerm = cplex.linearNumExpr();
        decTerm.addTerm(1, bidderTwoDecisions[i]);
        cplex.addLe(probTerm, decTerm);
    }       
    //expression (23)    
    for (int i = 0; i < numTypesOne; i++)
    {
        IloLinearNumExpr expr= cplex.linearNumExpr();
        expr.addTerm(1, typeUtilsOne[i]);
        for (int k = 0; k < numTypesOne; k++)
        {
            IloLinearNumExpr expr2 = cplex.linearNumExpr();
            expr.addTerm(1, bidderOneEUs[i * numTypesOne + k]);
            cplex.addLe(expr2, expr);
        }
    }
    for (int i = 0; i < numTypesTwo; i++)
    {
        IloLinearNumExpr expr= cplex.linearNumExpr();
        expr.addTerm(1, typeUtilsTwo[i]);
        for (int k = 0; k < numTypesTwo; k++)
        {
            IloLinearNumExpr expr2 = cplex.linearNumExpr();
            expr.addTerm(1, bidderTwoEUs[i * numTypesTwo + k]);
            cplex.addLe(expr2, expr);
        }
    }
    //expression (22)
    for (int i = 0; i < numTypesOne; i++)
    {
        for (int k = 0; k < numTypesOne; k++)
        {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1, bidderOneEUs[i * numTypesOne + k]);
            
            IloLinearNumExpr expr2 = cplex.linearNumExpr();
            expr2.addTerm(1, typeUtilsOne[i]);        
            expr2.addTerm(Double.MAX_VALUE/2, bidderOneDecisions[i * numTypesOne + k]);
            expr2.setConstant(expr2.getConstant() - Double.MAX_VALUE/2);
        
            cplex.addGe(expr, expr2);
        }
    }
    for (int i = 0; i < numTypesTwo; i++)
    {
        for (int k = 0; k < numTypesTwo; k++)
        {
            IloLinearNumExpr expr = cplex.linearNumExpr();
            expr.addTerm(1, bidderTwoEUs[i * numTypesTwo + k]);
            
            IloLinearNumExpr expr2 = cplex.linearNumExpr();
            expr2.addTerm(1, typeUtilsTwo[i]);        
            expr2.addTerm(Double.MAX_VALUE/2, bidderTwoDecisions[i * numTypesTwo + k]);
            expr2.setConstant(expr2.getConstant() - Double.MAX_VALUE/2);
        
            cplex.addGe(expr, expr2);
        }
    }
    //expression (20)
    for (int i = 0; i < numTypesOne; i++)
    {
        for (int k = 0; k < numTypesOne; k++)
        {
            IloLinearNumExpr eu = cplex.linearNumExpr();
            for (int q = 0; q < numTypesTwo; q++)
            {
                double probType = probType(2, q, dist);
                for (int w = 0; w < numTypesTwo; w++)
                {
                    double utility = solution.probAgentOne[k * numTypesTwo + w] * ((double)i) - solution.paymentsAgentOne[k * numTypesTwo + w];
                    eu.addTerm(utility * probType, bidderTwoProbs[q * numTypesTwo + w]);
                }
            }
            
            IloLinearNumExpr bidderOneEU = cplex.linearNumExpr();
            bidderOneEU.addTerm(1, bidderOneEUs[i * numTypesOne + k]);
            cplex.addEq(eu, bidderOneEU);
        }
    }/*
    //expression (21)
    for (int i = 0; i < numTypesTwo; i++)
    {
        for (int k = 0; k < numTypesTwo; k++)
        {
            IloLinearNumExpr eu = cplex.linearNumExpr();
            for (int q = 0; q < numTypesOne; q++)
            {
                double probType = probType(1, q, dist);
                for (int w = 0; w < numTypesOne; w++)
                {
                    double utility = solution.probAgentTwo[w * numTypesTwo + k] * ((double)i) - solution.paymentsAgentTwo[w * numTypesTwo + k];
                    eu.addTerm(utility * probType, bidderOneProbs[q * numTypesOne + w]);
                }
            }
            
            IloLinearNumExpr bidderTwoEU = cplex.linearNumExpr();
            bidderTwoEU.addTerm(1, bidderTwoEUs[i * numTypesTwo + k]);
            cplex.addEq(eu, bidderTwoEU);
        }
    }
   */

    //objective function
    IloLinearNumExpr obj = cplex.linearNumExpr();
    for (int i = 0; i < numTypesOne; i++)
    {
        double probType = probType(1, i, dist);
        for (int k = 0; k < numTypesOne; k++)
        {
            obj.addTerm(probType, bidderOneEUs[i * numTypesOne + k]);
        }
    }
    for (int i = 0; i < numTypesTwo; i++)
    {
        double probType = probType(2, i, dist);
        for (int k = 0; k < numTypesTwo; k++)
        {
            obj.addTerm(probType, bidderTwoEUs[i * numTypesTwo + k]);
        }
    }
    System.out.println("STUFF");
    //IloObjective objective = cplex.addMaximize(obj);
    boolean b = cplex.solve();
    System.out.println("SOLVED = " + b);
    double objVal = 5;//cplex.getObjValue();
    cplex.end();       
    //System.out.println(objVal);
    return objVal; 
  }

}
