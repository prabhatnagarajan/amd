import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.util.*;
import java.util.*;
import java.io.*;
//import java.awt.Pair;
public class BayesianStatistics
{
    public static class Sextuple
    {
        public ArrayList<Pair<Integer, Double>> lowerOne;
        public ArrayList<Pair<Integer, Double>> upperOne;
        public ArrayList<Pair<Integer, Double>> lowerTwo;
        public ArrayList<Pair<Integer, Double>> upperTwo;
        public ArrayList<Pair<Integer, Double>> avg;
        public List<Pair<Integer, Double>> actual;

        public Sextuple(ArrayList<Pair<Integer, Double>> lowerOne, ArrayList<Pair<Integer, Double>> upperOne, ArrayList<Pair<Integer, Double>> lowerTwo, ArrayList<Pair<Integer, Double>> upperTwo, ArrayList<Pair<Integer, Double>> avg, List<Pair<Integer, Double>> actual)
        {
            this.lowerOne = lowerOne;
            this.upperOne = upperOne;
            this.lowerTwo = lowerTwo;
            this.upperTwo = upperTwo;
            this.avg = avg;
            this.actual = actual;
        }
    }

    public Sextuple sextuple;
    public int numTypesOne;
    public int numTypesTwo;
    public Double[] agentOneTypes;
    public Double[] agentTwoTypes;
    public long learnTime;
    public BayesianStatistics(String typeFile, String distFile, int actualSamples, int numDSamples, double confidence) throws IOException
    {
        if (confidence > 1 || confidence < 0)
        {
            throw new IOException();
        }
        Scanner typeReader = new Scanner(new File(typeFile));
        String line = typeReader.nextLine();
        String line2 = typeReader.nextLine();
        agentOneTypes = getAgentOneTypes(line);
        agentTwoTypes = getAgentTwoTypes(line2);
        numTypesOne = agentOneTypes.length;
        numTypesTwo = agentTwoTypes.length;
        EnumeratedDistribution<Integer> typeDistribution = getDistribution(distFile);
        long before = System.currentTimeMillis();
        sextuple = learn2(typeDistribution, actualSamples, numDSamples, confidence);
        long after = System.currentTimeMillis();
        learnTime = after - before;
    }

    public Sextuple learn(EnumeratedDistribution<Integer> actual, int actualSamples, int numDSamples, double confidence)
    {
        /*
        STEPS:
        Set Dirichlet to alpha vector = 1
        Obtain Samples from the actual distribution
        Sample from the Dirichlet Distribution and obtain confidence interval, return lower and upper bound
        */
        /*for (Pair<Integer, Double> pair : actual.getPmf())
        {
            System.out.println(pair.getSecond());
        }
        */
        int length = actual.getPmf().size();
        double []alpha= new double[length];
        for (int i = 0; i < alpha.length; i++)
        {
            alpha[i] = 1;
        }
        DirichletDistribution dist = new DirichletDistribution(alpha);
        int numSamples = actualSamples;
        Integer []actualSample = new Integer[numSamples];
        actual.sample(numSamples, actualSample);
        double []observations = new double[length];
        for (int i = 0; i < actualSample.length; i++)
        {
            int res = actualSample[i];
            observations[res] = observations[res] + 1;
        }
        dist.updateHyperParameter(observations);
        //POTENTIALLY CHANGE, TESTING NEW DIRICHLET SAMPLING
        // double [][]dirichletSamples = dist.sample(numDSamples);
       // long before = System.currentTimeMillis();
       // double [][]dirichletSamples = dist.sample(numDSamples);
       // long after = System.currentTimeMillis();
        long before = System.currentTimeMillis();
        Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> lowup = dist.sample(numDSamples, confidence); 
        long after = System.currentTimeMillis();
        System.out.println("Dirichlet Sampling of " + numDSamples + " samples takes " + (after- before) + " Milliseconds.");
        //Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> lowup = confidenceInterval(confidence, dirichletSamples); 
        //HERE I SHOULD CREATE AN AGGREGATE LOWER AND UPPER BOUND FROM EACH INDIVIDUAL DIRICHLET 
        return new Sextuple(lowup.getFirst(), lowup.getSecond(), null, null, getMean(dist), actual.getPmf());
    }

    public Sextuple learn2(EnumeratedDistribution<Integer> actual, int actualSamples, int numDSamples, double confidence)
    {

        int length = actual.getPmf().size();
        int numSamples = actualSamples;
        Integer []actualSample = new Integer[numSamples];

        double []alphaG= new double[length];
        for (int i = 0; i < alphaG.length; i++)
        {
            alphaG[i] = 1;
        }
        DirichletDistribution dist = new DirichletDistribution(alphaG);
        actual.sample(numSamples, actualSample);
        double []observations = new double[length];
        for (int i = 0; i < actualSample.length; i++)
        {
            int res = actualSample[i];
            observations[res] = observations[res] + 1;
        }
        dist.updateHyperParameter(observations);

        DirichletDistribution []condDistsOne = new DirichletDistribution[numTypesTwo];
        //One Given Two
        for (int i = 0; i < condDistsOne.length; i++)
        {
            double []alpha= new double[numTypesOne];
            for (int k = 0; k < alpha.length; k++)
            {
                alpha[k] = 1;
            }
            condDistsOne[i] = new DirichletDistribution(alpha);
            double obs[] = new double[numTypesOne];
            for (int k = 0; k < alpha.length; k++)
            {
                obs[k] = observations[k * numTypesTwo + i];
            }
            condDistsOne[i].updateHyperParameter(obs);
        }

        DirichletDistribution []condDistsTwo = new DirichletDistribution[numTypesOne];
        //Two Given One
        for (int i = 0; i < condDistsTwo.length; i++)
        {
            double []alpha= new double[numTypesTwo];
            for (int k = 0; k < alpha.length; k++)
            {
                alpha[k] = 1;
            }
            condDistsTwo[i] = new DirichletDistribution(alpha);
            double obs[] = new double[numTypesTwo];
            for (int k = 0; k < alpha.length; k++)
            {
                obs[k] = observations[i * numTypesTwo + k];
            }
            condDistsTwo[i].updateHyperParameter(obs);
        }

        ArrayList<Pair<Integer, Double>> oneLower = new ArrayList<Pair<Integer,Double>>(); 
        ArrayList<Pair<Integer, Double>> oneUpper = new ArrayList<Pair<Integer,Double>>();
        ArrayList<Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>>> aggregateOne= new ArrayList<Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>>>();
        for (int i = 0; i < condDistsOne.length; i++)
        {
            Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> lowupOne = condDistsOne[i].sample(numDSamples, confidence);
            aggregateOne.add(lowupOne);
        }

        for (int i = 0; i < numTypesOne; i++)
        {
            for (int k = 0; k < numTypesTwo; k++)
            {
                oneLower.add(aggregateOne.get(k).getFirst().get(i));
                oneUpper.add(aggregateOne.get(k).getSecond().get(i));
            }
        }

        ArrayList<Pair<Integer, Double>> twoLower = new ArrayList<Pair<Integer,Double>>(); 
        ArrayList<Pair<Integer, Double>> twoUpper = new ArrayList<Pair<Integer,Double>>();
        ArrayList<Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>>> aggregateTwo= new ArrayList<Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>>>();
        for (int i = 0; i < condDistsTwo.length; i++)
        {
            Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> lowupTwo = condDistsTwo[i].sample(numDSamples, confidence);
            aggregateTwo.add(lowupTwo);
        }

        for (int i = 0; i < numTypesOne; i++)
        {
            for (int k = 0; k < numTypesTwo; k++)
            {
                twoLower.add(aggregateTwo.get(i).getFirst().get(k));
                twoUpper.add(aggregateTwo.get(i).getSecond().get(k));
            }
        }

        Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> lowup = dist.sample(numDSamples, confidence); 
        //HERE I SHOULD CREATE AN AGGREGATE LOWER AND UPPER BOUND FROM EACH INDIVIDUAL DIRICHLET 
        return new Sextuple(oneLower, oneUpper, twoLower, twoUpper, getMean(dist), actual.getPmf());
    }
    
    //confidence level should be between 0 and 1, e.g. 0.95
    public Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> confidenceInterval(double confidenceLevel, double [][]dirichletSamples)
    {
        double significance = 1 - confidenceLevel;
        int length = dirichletSamples[0].length;
        double [][]rotatedMatrix = new double[length][dirichletSamples.length];
        for (int i = 0; i < dirichletSamples.length; i++)
        {
            for (int k = 0; k < dirichletSamples[i].length; k++)
            {
                rotatedMatrix[k][i] = dirichletSamples[i][k];
            }
        }
        double []lower = new double[length];
        double []upper = new double[length];
        for (int i = 0; i < rotatedMatrix.length; i++)
        {
            double []varDist = rotatedMatrix[i];
            Arrays.sort(varDist);
            int lowIndex = (int) ((significance/2) * (varDist.length - 1));
            int highIndex = (int) ((1 - significance/2) * (varDist.length - 1));
            lower[i] = varDist[lowIndex];
            upper[i] = varDist[highIndex];
        }
        ArrayList<Pair<Integer, Double>> lowerBound = new ArrayList<Pair<Integer, Double>>();
        ArrayList<Pair<Integer, Double>> upperBound = new ArrayList<Pair<Integer, Double>>();
        for (int i = 0; i < lower.length; i++)
        {
            lowerBound.add(new Pair<Integer, Double>(i, lower[i]));
            upperBound.add(new Pair<Integer, Double>(i, upper[i]));
        }

        return new Pair<ArrayList<Pair<Integer, Double>>, ArrayList<Pair<Integer, Double>>> (lowerBound, upperBound);
    }

    public ArrayList<Pair<Integer, Double>> getMean(DirichletDistribution dirichlet)
    {
        double sum = 0;
        double []mean = new double[dirichlet.alpha.length];
        for (int i = 0; i < dirichlet.alpha.length; i++)
        {
            sum += dirichlet.alpha[i]; 
        }

        for (int i = 0; i < mean.length; i++)
        {
            mean[i] = dirichlet.alpha[i]/sum;
        }
        
        ArrayList<Pair<Integer, Double>> distribution = new ArrayList<Pair<Integer, Double>>();
        for (int i = 0; i < mean.length; i++)
        {
            distribution.add(new Pair<Integer, Double>(i, mean[i]));
        }
        return distribution;
    }

	public Double[] getAgentOneTypes(String line)
	{
		Scanner type1 = new Scanner(line);
		ArrayList<Double> list1 = new ArrayList<Double>();
		int count1 = 0;
		while (type1.hasNext())
		{
			list1.add(type1.nextDouble());
			count1++;
		}
		Double agentOneTypes[] = new Double[count1];
		agentOneTypes = list1.toArray(agentOneTypes);
		return agentOneTypes;
	}

	public Double[] getAgentTwoTypes(String line2)
	{
		Scanner type2 =  new Scanner(line2);
		ArrayList<Double> list2 = new ArrayList<Double>();
		int count2 = 0;
		while (type2.hasNext())
		{
			list2.add(type2.nextDouble());
			count2++;
		}
		Double agentTwoTypes[] = new Double[count2];
		agentTwoTypes = list2.toArray(agentTwoTypes);
		return agentTwoTypes;
	}

	public EnumeratedDistribution<Integer> getDistribution(String filename) throws IOException
	{
		Scanner input = new Scanner(new File(filename));
		ArrayList<Pair<Integer, Double>> list = new ArrayList<Pair<Integer, Double>>();
		int count = 0;
		while (input.hasNext())
		{
			list.add(new Pair<Integer, Double>(count, Math.max(input.nextDouble(), 0)));
			count++;
		}
		return new EnumeratedDistribution<Integer>(list);
	}
}
