// See https://aka.ms/new-console-template for more information


using System; 
using System.Collections.Generic;
using System.Text;
using Npgsql;

namespace MyApp
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Sum Twelve method:");
            Console.WriteLine(SumTwelve());
            Console.WriteLine("Box Muller method:");
            Console.WriteLine(BoxMuller()[0]);
            Console.WriteLine("Polar Rejection method:");
            Console.WriteLine(PolarRejection()[0]);
            Console.WriteLine("Here is the pair of correlated numbers.");
            double[] correlatedPair = Correlator(BoxMuller()[0],PolarRejection()[0],.23);
            Console.WriteLine(correlatedPair[0].ToString() + " and " + correlatedPair[1].ToString());

        }

        static double SumTwelve( )
        {
            // makes a draw from the normal distribution according to the sum twelve method.
            Random r = new Random();
            List<double> uniformDraws = new List<double>();
            for(int i = 0; i < 12; i++ )
                {
                    uniformDraws.Add(r.NextDouble());
                }
            double normalDraw = uniformDraws.Sum() - 6.0;
            return normalDraw;
        }
        
        static double[] BoxMuller()
        {
            // makes a draw from the normal distribution via the box muller method
            // RETURNS A PAIR OF NUMBERS

            Random r = new Random();
            double x1 = r.NextDouble();
            double x2 = r.NextDouble();
            double z1 = Math.Sqrt(-2*Math.Log(x1))*Math.Cos(2*Math.PI*x2);
            double z2 = Math.Sqrt(-2*Math.Log(x1))*Math.Sin(2*Math.PI*x2);
            double[] normalDraws = {z1,z2};
            return normalDraws;
        }

        static double[] PolarRejection()
        {   
            // makes a draw from the normal distribution according to the polar rejection method
            // RETURNS A PAIR OF NUMBERS
            Random r = new Random();
            double w = 2;
            double x1 = r.NextDouble();
            double x2 = r.NextDouble();
            while(w > 1)
            {
                x1 = r.NextDouble();
                x2 = r.NextDouble();
                w = Math.Pow(x1,2) + Math.Pow(x2,2);
            }
            double c = Math.Sqrt((-2*Math.Log(w))/w);
            double z1 = c*x1;
            double z2 = c*x2;
            double[] normalDraws = {z1,z2};
            return normalDraws;

        }

        static double[] Correlator(double epsilon1, double epsilon2, double rho)
        {
            // correlates two uncorrelated draws from the normal distribution.
            // RETURNS A PAIR OF NUMBERS
            double z2 = rho*epsilon1 + Math.Sqrt(1 - Math.Pow(rho,2))*epsilon2;
            double[] correlatedPair = {epsilon1,z2};
            return correlatedPair;
        }


    }
}

