import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.util.*;
import java.util.*;
import java.io.*;
public class TwoBidderDistributions
{

	public double distribution[][];
	public double oneGivenTwo[][];
	public double twoGivenOne[][];
	public EnumeratedDistribution<Integer> oneGivenTwoDist;
	public EnumeratedDistribution<Integer> twoGivenOneDist;
	public int numTypes1;
	public int numTypes2;

	public TwoBidderDistributions(double matrix[][], int numTypes1, int numTypes2)
	{
		this.distribution = matrix;
		this.numTypes1 = numTypes1;
		this.numTypes2 = numTypes2;
		obtainOneGivenTwo();
		obtainTwoGivenOne();
		getEnumeratedDistributionOne();
		getEnumeratedDistributionTwo();
	}

	private void obtainOneGivenTwo()
	{
		oneGivenTwo = new double[numTypes1][numTypes2];
		for (int i = 0; i < numTypes2; i++)
		{
			double sum = sumColumn(distribution, i);
			for (int k = 0; k < numTypes1; k++)
			{
				oneGivenTwo[k][i] = distribution[k][i]/sum;
			}
		}
	}

	private void obtainTwoGivenOne()
	{
		twoGivenOne = new double[numTypes1][numTypes2];
		for (int i = 0; i < numTypes1; i++)
		{
			double sum = sumRow(distribution, i);
			for (int k = 0; k < numTypes2; k++)
			{
				twoGivenOne[i][k] = distribution[i][k]/sum;
			}
		}	
	}

	private double sumColumn(double matrix[][], int columnNum)
	{
		double sum = 0;
		for (int i = 0; i < matrix.length; i++)
		{
			sum += matrix[i][columnNum];
		}
		return sum;
	}

	private double sumRow(double matrix[][], int rowNum)
	{
		double sum = 0;
		for (int i = 0; i < matrix[rowNum].length; i++)
		{
			sum += matrix[rowNum][i];
		}
		return sum;
	}

	private void getEnumeratedDistributionOne()
	{
		ArrayList<Pair<Integer, Double>> list = new ArrayList<Pair<Integer, Double>>();
		for (int i = 0; i < numTypes1; i++)
		{
			for (int k = 0; k < numTypes2; k++)
			{
				list.add(new Pair<Integer, Double>(i * numTypes2 + k, Math.max(oneGivenTwo[i][k], 0)));
			}
		}
		oneGivenTwoDist = new EnumeratedDistribution<Integer>(list);
	}

	private void getEnumeratedDistributionTwo()
	{
		ArrayList<Pair<Integer, Double>> list = new ArrayList<Pair<Integer, Double>>();
		for (int i = 0; i < numTypes1; i++)
		{
			for (int k = 0; k < numTypes2; k++)
			{
				list.add(new Pair<Integer, Double>(i * numTypes2 + k, Math.max(twoGivenOne[i][k], 0)));
			}
		}
		twoGivenOneDist = new EnumeratedDistribution<Integer>(list);
	}
}