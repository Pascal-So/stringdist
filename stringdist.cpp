#include<bits/stdc++.h>
using namespace std;

bool trace = false;
bool verbose = true;

// implementing an algorithm to map out
// words on a plane, grouping similar
// ones together (similar as defined
// by low levenshtein distance).


// levenshtein string distance
int stringdist(string stringA, string stringB)
{

  int lengthA = stringA.size();
  int lengthB = stringB.size();

  // only store 2 rows and swap them on every step
  vector<vector<int> > rows (2, vector<int> (lengthA+1, 0));

  for(int i = 0; i < lengthA; ++i){
    rows[0][i] = i;
  }
  
  for(int j = 1; j <= lengthB; ++j){
    rows[1][0]=j;

    for(int i = 1; i <= lengthA; ++i){

      // min of insertStep and deleteStep
      rows[1][i] = min(rows[0][i]+1,
		       rows[1][i-1]+1
		       );

      // if chars are equal, the cost to replace
      // one by the other is 0, otherwise it's 1
      int replaceCost = stringA[i-1]==stringB[j-1]? 0 : 1;

      rows[1][i] = min(rows[1][i],
		       replaceCost + rows[0][i-1]
		       );
    }
    
    
    swap(rows[0], rows[1]);
  }

  return rows[0][lengthA];
}

vector<vector<int> > buildDistanceMatrix(int n, vector<string> strings)
{

  vector<vector<int> > dists (n, vector<int> (n, 0));
  
  // O(n^2) populate distance matrix
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < i; ++j){
      int dist = stringdist(strings[i], strings[j]);
      dists[i][j] = dist;
      dists[j][i] = dist;
    }

    // print out the current percentage at 20% steps
    if (verbose && n >= 10){
      if(i % (n/5) == 0 && i >= n/10){
	cerr<<i*100/n << "% done \n";
      }
    }
  }

  return dists;
}

double drand(){
  return (double)rand() / RAND_MAX;
}

void randomizeEmbedding(int n, vector<pair<double, double> > & embedding)
{
  for(int i=0; i < n; ++i){
    embedding[i].first = drand();
    embedding[i].second = drand();
  }
}



double dist(pair<double, double> pA, pair<double, double> pB){
  double dx = pB.first - pA.first;
  double dy = pB.second - pA.second;

  return sqrt(dx*dx + dy*dy);
}

pair<double, double> costDerivative( int n,
				     vector<pair<double, double> > & embedding,
				     vector<vector<int> > & distanceMatrix,
				     int point,
				     double distScale)
{

  // derivative of sammon mapping
  
  // derivatives will be summed up on these variables,
  // considering every other node individually
  double dx = 0;
  double dy = 0;
  
  for(int i = 0; i < n; ++i){
    if(i==point){
      continue;
    }


    // calculate the magnitude of the derivative first
    double distOrig = (double) distanceMatrix[i][point] * distScale;
    double distMapping = dist(embedding[i], embedding[point]);


    // f  = (dO - dM)^2 / dO  the cost function
    // f' = 2*dM / dO - 2     derivative of the cost function
    double derivMag = 2*distMapping / distOrig - 2;

    // split the derivative in to its components
    // and add it to the sum of contributors 
    dx += derivMag * (embedding[i].first - embedding[point].first) / distOrig;
    dy += derivMag * (embedding[i].second - embedding[point].second) / distOrig;
  }

  return {dx, dy};
}


// modifies values in argument 'embedding'
void descentPlaneEmbedding(int n,
			   vector<pair<double, double> > & embedding,
			   vector<vector<int> > & distanceMatrix,
			   int iterations,
			   double stepFactor,
			   double distScale
			   )
{
  
  for(int it = 0; it < iterations; ++it){

    vector<pair<double, double> > offsets (n);
    
    for(int i = 0; i < n; ++i){
      pair<double, double> derivative = costDerivative(n, embedding, distanceMatrix, i, distScale);

      derivative.first *= stepFactor;
      derivative.second *= stepFactor;

      // store derivative and apply it after we've
      // calculated all derivatives, otherwise we'd
      // be considering different positions for
      // the later calculations.
      offsets[i] = derivative;

    }

    for(int i = 0; i < n; ++i){
      // apply all offsets
      embedding[i].first += offsets[i].first;
      embedding[i].second += offsets[i].second;
    }

    // this could be turned on by setting the `trace`
    // variable at the top to true. This will print out
    // not just the final embedding, but the positions
    // of the points at every timestep, useful for
    // showing trails of the motion in the final plot
    if(trace){
      for(auto & e : embedding)
	cout<< e.first << " " << e.second<<"\n"; 
      cout<<"\n\n";
    }

    // print out the current percentage at 10% steps
    if (verbose && iterations >= 10){
      if(it % (iterations/10) == 0 && it >= iterations/10){
	cerr<<it*100/iterations << "% done \n";
      }
    }
  }
  
}

int main(int argc, char* argv[])
{

  int iterations;
  if(argc > 1){
    istringstream(argv[1]) >> iterations;
  }else{
    iterations = 500;
  }
  
  double stepFactor = 0.0002;
  double distScale = 0.1;
  
  
  srand(time(NULL));
  
  vector<string> strings;
  string s;
  while(getline(cin, s)){
    strings.push_back(s);
  }

  int n = strings.size();
  
  if(verbose) cerr<<"Building " << n << "*" << n << " distance matrix..\n";
  
  vector<vector<int> > distanceMatrix = buildDistanceMatrix(n, strings);
  
  if(verbose) cerr<<"Distance matrix done, starting embedding with " << iterations << " iterations.\n";
  
  vector<pair<double, double> > embedding (n);

  // initialize the embedding with random values
  randomizeEmbedding(n, embedding);

  // adjust the embedding iteratively by gradient descent
  descentPlaneEmbedding(n, embedding, distanceMatrix, iterations, stepFactor, distScale);


  // print out the final embedding as "x y"
  for(auto & e : embedding)
    cout<< e.first << " " << e.second<<"\n"; 

}
