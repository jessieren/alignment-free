
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include<cmath>
#include <vector>
#include <sstream>
#include <string>
#include <stdlib.h> 
//strtol(s.c_str(),0,10);
using namespace std;

# include <stdio.h>
# include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include<algorithm>
#include <tr1/unordered_map>
using namespace std::tr1;
using namespace std;

int ZI = 4;

struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};

struct KmerInfo
{
	//receptacle of kmer counts
	unordered_map<unsigned long, unsigned long> HashTable;

	unordered_map<unsigned long, unsigned long> HashTableK_1;
	// kmer length = k - 2
	unordered_map<unsigned long, unsigned long> HashTableK_2;

	unordered_map<unsigned long, unsigned long> HashTableOrder;
	// kmer length = order + 1
	unordered_map<unsigned long, unsigned long> HashTableOrder_1;

	//receptacle of Pw(probability of a kmer word) 
	unordered_map<unsigned long, SCIENTIFIC_NUMBER > HashPw;

	unsigned long totalKmer;
	unsigned long totalK_1;
	unsigned long totalK_2;
	unsigned long totalOrder;
	unsigned long totalOrder_1;
};


// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
  return dSci.value * pow(10, dSci.factor);
}


// TransToScientific : trans a read number to a SCIENTIFIC NUMBER
SCIENTIFIC_NUMBER TransToScientific(double dReal)
{
	SCIENTIFIC_NUMBER sciTemp;
	int count;
  
	if( dReal==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
	}
	else if(dReal>10.0 || dReal<-10.0)
	{
		count=0;
		while(dReal>10.0 || dReal<-10.0)
		{
			dReal /=10.0;
			count++;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else if( dReal<1.0 && dReal>-1.0)
	{
		count=0;
		while( dReal<1.0 && dReal>-1.0 )
		{
			dReal *=10.0;
			count--;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else
	{
		sciTemp.value=dReal;
		sciTemp.factor=0;
	}
  
	return sciTemp;
}


SCIENTIFIC_NUMBER SciNegative(SCIENTIFIC_NUMBER Sci)
{
  SCIENTIFIC_NUMBER sciNeg;
  sciNeg.value = - Sci.value; sciNeg.factor = Sci.factor;
  return sciNeg;
}

SCIENTIFIC_NUMBER SciInverse(SCIENTIFIC_NUMBER Sci)
{
  SCIENTIFIC_NUMBER sciInv;
  sciInv.value = 1/Sci.value; sciInv.factor = -Sci.factor;
  return sciInv;
}

// SciMultiple : Multiplication of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
  //    cout << "SciMultiple " << endl;
  
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int count;
  
	if( left.value==0.0 || right.value==0.0 )
	{
    //        cout << "Both 0 " << endl;
    
		sciTemp.value=0.0;
		sciTemp.factor=0;
    
		return sciTemp;
	}
  
	// now both left and right element are nonzero
	dTemp=left.value * right.value;
  
  //    cout << "left.value " << left.value << endl;
  //    cout << "right.value " << right.value << endl;
  //    cout << "dTemp " << dTemp << endl;
  
  
	if( dTemp>10.0 || dTemp<-10.0 )
	{
    
    //        cout << "10 < dTemp or dTemp < -10 " << endl;
    
		count=0;
		while(dTemp>10.0 || dTemp<-10.0 )
		{
			dTemp /=10.0;
			count++;
		}
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
	}
	else if( dTemp<1.0 && dTemp>-1.0)
	{
    //        cout << "dTemp < 1 or dTemp > -1 " << dTemp << endl;
    
		count=0;
		while( dTemp<1.0 && dTemp>-1.0 )
		{
			dTemp *=10.0;
			count--;
		}
    
    //        cout << "count " << count << endl;
    
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
    
    //        cout << "sciTemp " << sciTemp.value << " " << sciTemp.factor << endl;
	}
	else
	{
    
    //        cout << "dTemp normal " << dTemp << endl;
    
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor;
	}
  
	return sciTemp;
}


// SciMultiple : Multiplication between a SCIENTIFIC NUMBERS and a read number
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
  
	sciTemp=TransToScientific(right);
	sciTemp=SciMultiple(left,sciTemp);
  
	return sciTemp;
}





// SciAddition : addition of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int i,count;
  
	if( left.value==0.0 || right.value==0.0 )
	{
//    cout << "this step"  << endl;
		if( left.value==0.0 )
    {
			return right;
    }else{
			return left;
    }
	}
  
	// now the two element are both non zero
	if( left.factor>=right.factor)
	{
		// left element is larger than right element
		dTemp=right.value;
		for(i=0;i<(left.factor-right.factor);i++)
			dTemp /=10.0;
		dTemp +=left.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}
    
		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while(dTemp>10.0 || dTemp<-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor;
		}
		return sciTemp;
	}
	else
	{
		// right element  is larger than left element
		dTemp=left.value;
		for(i=0;i<(right.factor-left.factor);i++)
			dTemp /=10.0;
		dTemp +=right.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}
    
		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while( dTemp>10.0 || dTemp <-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor;
		}
		return sciTemp;
	}
}


// SciAddition : addition between a SCIENTIFIC NUMBERS and a real number
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
  
	sciTemp=TransToScientific(right);
	sciTemp=SciAddition(left,sciTemp);
  
	return sciTemp;
}

// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciPow(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	double dTemp;
	int iTemp;
  /*
   if(left.value==0.0 )
   {
   printf("the base of the power is nagative\n");
   exit(1);
   }
   */
	if(left.value==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
		return sciTemp;
	}
  
	dTemp=(log10(fabs(left.value))+left.factor)*right;
  
	if( dTemp>0.0 )
		iTemp=int(ceil(dTemp));    //ceil(a)是求不小于a的最小整数。floor(a)表示求不大于a的最大整数
	else
		iTemp=int(floor(dTemp));
	sciTemp.factor=iTemp;
	sciTemp.value=pow(10.0,dTemp-iTemp);
  
	return sciTemp;
}





// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciLn(SCIENTIFIC_NUMBER sci)
{
  SCIENTIFIC_NUMBER sciTemp = TransToScientific( log( sci.value ) + sci.factor * log(10) ) ;
  
  return sciTemp;
}





void printFour(vector<int> four)
{
  cout << "print ";
  for(int it=0; it<four.size(); it++)
  {
    cout << four[it] << "," ;
  }
  cout << endl;

}



vector<int> ten2four(unsigned long ten, int k)
{
  vector<int> four (k,0);
  unsigned long tmp = ten; 
  int currentPos = k-1;
  for(int currentPos = k-1; currentPos >=0; --currentPos)
  {
    four[currentPos]=tmp%ZI;
    tmp/=ZI; 
  }
  //while(tmp>=(ZI-1)) {four[currentPos]=tmp%ZI;tmp/=ZI;currentPos--; }
  //four[currentPos] = tmp;
  return four;
}


int four2ten(vector<int> four, int k)
{
  unsigned long ten = 0;
  for(int currentPos=(k-1); currentPos >= 0; --currentPos)
  { 
    int tmp = four[currentPos] * pow(ZI,(k-1 - currentPos));
    ten = ten + tmp;
    //cout << currentPos << " " << ten << endl;
  }
  return ten;
  
}


void loadKmerCountHash(string currentKmerFilePathName, KmerInfo* speciesKmerInfo, string kmerUsage)
{
	ifstream fin(currentKmerFilePathName.c_str());
  string currentKmerLine;
  while(getline(fin, currentKmerLine, '\n'))
  { 
    unsigned long currentKmerID; unsigned long currentKmerCount;
    std::istringstream ss(currentKmerLine);
    std::string token;
    int colCount = 0;
    while(std::getline(ss, token, ',')) {
      colCount++;
      if(colCount == 1){
        currentKmerID = atoi(token.c_str());
      }else{
        currentKmerCount = atoi(token.c_str());
      }//std::cout << token << '\n';
    }
    //cout << "kmerID:" << currentKmerID << " kmerCount:" << currentKmerCount << endl;
    if(kmerUsage == "kmerCount"){
		speciesKmerInfo->HashTable[currentKmerID] = currentKmerCount;
		speciesKmerInfo->totalKmer += currentKmerCount;
    }else if(kmerUsage == "kmerCount_1"){
		speciesKmerInfo->HashTableK_1[currentKmerID] = currentKmerCount;
		speciesKmerInfo->totalK_1 += currentKmerCount;
    }else if(kmerUsage == "kmerCount_2"){
		speciesKmerInfo->HashTableK_2[currentKmerID] = currentKmerCount;
		speciesKmerInfo->totalK_2 += currentKmerCount;
    }else if(kmerUsage == "kmerOrder"){
		speciesKmerInfo->HashTableOrder[currentKmerID] = currentKmerCount;
		speciesKmerInfo->totalOrder += currentKmerCount;
    }else if(kmerUsage == "kmerOrder+1"){
		speciesKmerInfo->HashTableOrder_1[currentKmerID] = currentKmerCount;
		speciesKmerInfo->totalOrder_1 += currentKmerCount;
    }else{cout << "ERROR: wrong kmerUsage" << endl;}
    
  }
  return;
  
}



vector<int> reverseFour(vector<int> Four)
{
  vector<int> reverseFour(Four.size(), 4);
  for(int revPos = 0; revPos < Four.size(); revPos++)
  { 
    reverseFour[revPos] = 3 - Four[Four.size()- 1 - revPos];
  }
  return reverseFour;
  
}
  


void pwIID(int ZI, int k, KmerInfo* speciesKmerInfo)
{
  // compute pw for each kmer word
  unordered_map<unsigned long,SCIENTIFIC_NUMBER > probIID;
  for(int index = 0; index < ZI; index++)
  { 


	  probIID[index] = TransToScientific(speciesKmerInfo->HashTableOrder_1[index] / double(speciesKmerInfo->totalOrder_1));
  }
  
  for(unsigned long ten = 0; ten < pow(ZI, k); ten++ )
  {
    SCIENTIFIC_NUMBER pw;
    pw.value = 1; pw.factor = 0;
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(ten, k);
    for(int pos = 0; pos < k; pos ++)
    {
      pw = SciMultiple(pw, probIID[currentKmerFour[pos]]);
    }
    
	speciesKmerInfo->HashPw[ten] = pw;
  }
}



void pwMC(int ZI, int k, int order, KmerInfo* speciesKmerInfo, vector<SCIENTIFIC_NUMBER>& iniProb, vector< vector<SCIENTIFIC_NUMBER> > &transMatrix)
{
  // compute the transition probability matrix
  int transColSize = ZI;
  int transRowSize = pow(ZI, order);
  //vector< vector<SCIENTIFIC_NUMBER> > transMatrix(transRowSize, vector<SCIENTIFIC_NUMBER>(transColSize));
  //vector<SCIENTIFIC_NUMBER> iniProb;
  for(int currentRow = 0; currentRow < transRowSize; currentRow++)
  {
    vector<int> currentRowFour = ten2four(currentRow, order);
	unsigned long countBelow = speciesKmerInfo->HashTableOrder[currentRow];
	double probBelow = countBelow / double(speciesKmerInfo->totalOrder);
    SCIENTIFIC_NUMBER probBelowSci = TransToScientific(probBelow);
    iniProb.push_back(probBelowSci);
    // %% it could be 0 ! %%
    // 20140910: the denominator cannot be 0
    SCIENTIFIC_NUMBER inv_probBelowSci;
    inv_probBelowSci.value = 0; inv_probBelowSci.factor = 0;
    if(probBelowSci.value != 0)
    {
      inv_probBelowSci = SciInverse(probBelowSci);
    }
    
    SCIENTIFIC_NUMBER rowSum;
    rowSum.value = 0; rowSum.factor = 0;
  
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      //cout << "row:" << currentRow << ", col:" << currentCol << endl;
      vector<int> currentRowColFour(currentRowFour);
      currentRowColFour.push_back(currentCol); 
      for(int it=0; it< currentRowColFour.size(); it++)
      {
        //cout << currentRowColFour[it];
      }
      //cout << endl;
      int currentRowColTen = four2ten(currentRowColFour, (order+1));
	  unsigned long countAbove = speciesKmerInfo->HashTableOrder_1[currentRowColTen];
	  double probAbove = countAbove / double(speciesKmerInfo->totalOrder_1);
      SCIENTIFIC_NUMBER probAboveSci = TransToScientific(probAbove);

      if( probBelow != 0 )
      {
        transMatrix[currentRow][currentCol] = SciMultiple(probAboveSci,inv_probBelowSci);
      }else{
        
        transMatrix[currentRow][currentCol] = TransToScientific(0);
      }
      
      rowSum = SciAddition(rowSum, transMatrix[currentRow][currentCol]);    
   }
  
    
    
    // normalize trans matrix
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      transMatrix[currentRow][currentCol] = SciMultiple(transMatrix[currentRow][currentCol], SciInverse(rowSum));
      //cout << TransToReal(transMatrix[currentRow][currentCol]) << endl;
      
    }
  }
  
  
  // compute pw for each kmer word
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++ )
  {
    SCIENTIFIC_NUMBER pw;
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << endl << endl;
    //printFour(currentKmerFour);
    vector<int> iniWordFour;
    for(int iniWordPos = 0; iniWordPos < order; iniWordPos ++)
    {
      iniWordFour.push_back(currentKmerFour[iniWordPos]);
    }
    unsigned long iniWordTen = four2ten(iniWordFour, iniWordFour.size());
    pw = iniProb[iniWordTen];
    //cout << "wordten: " << currentKmerTen << endl << "wordfour: ";
    //cout << "iniProb " << TransToReal(pw);
    
    for(int currentTransToPos = order; currentTransToPos < k; currentTransToPos++)
    {
      vector<int> currentTransFromWordFour;
      for(int transWordPos = currentTransToPos - order; transWordPos < currentTransToPos; transWordPos ++)
      {
        //cout << currentTransToPos << ", " << transWordPos << endl;
        //cout << currentKmerFour[transWordPos] << ",";
        currentTransFromWordFour.push_back(currentKmerFour[transWordPos]);
      }
      //printFour(currentTransFromWordFour);
      unsigned long currentTransFromWordTen = four2ten(currentTransFromWordFour, currentTransFromWordFour.size());
      pw = SciMultiple(transMatrix[currentTransFromWordTen][currentKmerFour[currentTransToPos]], pw);
      //cout << " trans from ";
      //printFour(currentTransFromWordFour);
      //cout << "to " << currentKmerFour[currentTransToPos] << " with prob " << TransToReal(transMatrix[currentTransFromWordTen][currentKmerFour[currentTransToPos]]) << endl;
      //cout << pw.value << ", " << pw.factor << endl;
    }
    
	speciesKmerInfo->HashPw[currentKmerTen] = pw;
  }
  
}


vector<double> D2C2computeNGS(int ZI, int k, KmerInfo* speciesKmerInfoA, KmerInfo* speciesKmerInfoB)
{
  // compute D2 statistics
  SCIENTIFIC_NUMBER D2, D2star, D2shepp;
  D2.value = 0; D2.factor = 0;  
  D2star.value = 0; D2star.factor = 0;
  D2shepp.value = 0; D2shepp.factor = 0;
  SCIENTIFIC_NUMBER C2_below[2];
  C2_below[0].value = 0; C2_below[0].factor = 0;
  C2_below[1].value = 0; C2_below[1].factor = 0;
  SCIENTIFIC_NUMBER C2star_below[2];
  C2star_below[0].value = 0; C2star_below[0].factor = 0;
  C2star_below[1].value = 0; C2star_below[1].factor = 0;
  SCIENTIFIC_NUMBER C2shepp_below[2];
  C2shepp_below[0].value = 0; C2shepp_below[0].factor = 0;
  C2shepp_below[1].value = 0; C2shepp_below[1].factor = 0;
  
  //cout << TransToReal(D2) << endl;
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
    
    unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
    //cout << currentKmerRevTen << endl;
    //printFour(currentKmerFour);
    //printFour(currentKmerFour);
    //printFour(currentKmerRevFour);
    
    SCIENTIFIC_NUMBER p_w[2], EX_w[2], X_w[2], X_w_tilde[2], X_w_tilde_var[2];
    
    for(int speciesID=0; speciesID<2; speciesID++)
    { 
		KmerInfo* speciesKmerInfo;
		if (speciesID == 0)
			speciesKmerInfo = speciesKmerInfoA;
		else
			speciesKmerInfo = speciesKmerInfoB;


		p_w[speciesID] = SciAddition(speciesKmerInfo->HashPw[currentKmerTen], speciesKmerInfo->HashPw[currentKmerRevTen]);
		EX_w[speciesID] = SciMultiple(p_w[speciesID], speciesKmerInfo->totalKmer);
		X_w[speciesID] = TransToScientific(speciesKmerInfo->HashTable[currentKmerTen] + speciesKmerInfo->HashTable[currentKmerRevTen]);
        X_w_tilde[speciesID] = SciAddition(X_w[speciesID], SciNegative(EX_w[speciesID]));
      
      // prepare for C2 compute
      C2_below[speciesID] = SciAddition(C2_below[speciesID], SciPow(X_w[speciesID], 2));
      // prepare for C2star compute
      // 20140910: the denominator cannot be 0
      if(EX_w[speciesID].value != 0)
      {
        X_w_tilde_var[speciesID] = SciMultiple(X_w_tilde[speciesID], SciInverse(SciPow(EX_w[speciesID],0.5)));
        C2star_below[speciesID] = SciAddition(C2star_below[speciesID], SciPow(X_w_tilde_var[speciesID], 2));
      }
      
      //cout << TransToReal(X_w[speciesID]) << endl;
    }
    
    //cout << "D2:" << D2.value << ", " << D2.factor << endl;
    //if(D2.value == 0.0){cout << "D2.value = 0.0" << endl;}
    D2 = SciAddition(SciMultiple(X_w[0], X_w[1]), D2);
    //cout << "X_w[0]:" << TransToReal(X_w[1]) << " X_w[1]:" << TransToReal(X_w[1]) << " D2:" << TransToReal(D2) << endl;
    
    SCIENTIFIC_NUMBER D2star_above = SciMultiple(X_w_tilde[0], X_w_tilde[1]);
    SCIENTIFIC_NUMBER D2star_below = SciPow(SciMultiple(EX_w[0], EX_w[1]), 0.5);
    if(D2star_below.value != 0)
    {
      D2star = SciAddition(D2star, SciMultiple(D2star_above, SciInverse(D2star_below)));
    }
    
    SCIENTIFIC_NUMBER D2shepp_below = SciPow(SciAddition(SciPow(X_w_tilde[0],2), SciPow(X_w_tilde[1],2)),0.5);
    // 20140910: the denominator cannot be 0
    if(D2shepp_below.value != 0)
    {
      D2shepp = SciAddition(D2shepp, SciMultiple(D2star_above, SciInverse(D2shepp_below)));
      // prepare for c2shepp compute
      C2shepp_below[0] = SciAddition(C2shepp_below[0], SciMultiple(SciPow(X_w_tilde[0],2), SciInverse(D2shepp_below)));
      C2shepp_below[1] = SciAddition(C2shepp_below[1], SciMultiple(SciPow(X_w_tilde[1],2), SciInverse(D2shepp_below)));
    }
    
  }
  
  SCIENTIFIC_NUMBER C2 = SciMultiple(D2, SciInverse(SciMultiple(SciPow(C2_below[0],0.5),SciPow(C2_below[1],0.5))));
  
  SCIENTIFIC_NUMBER C2star = SciMultiple(D2star, SciInverse( SciMultiple(SciPow(C2star_below[0], 0.5), SciPow(C2star_below[1], 0.5)) ));
  
  SCIENTIFIC_NUMBER C2shepp = SciMultiple(D2shepp, SciInverse( SciMultiple(SciPow(C2shepp_below[0],0.5),SciPow(C2shepp_below[1],0.5) )) );
  //cout << "d2shepp_below " << TransToReal(SciMultiple(SciPow(C2shepp_below[0],0.5),SciPow(C2shepp_below[1],0.5) )) << endl;
  

  vector<double> D2C2values;
  D2C2values.push_back(TransToReal(D2));
  D2C2values.push_back(TransToReal(C2));
  D2C2values.push_back(TransToReal(D2star));
  D2C2values.push_back(TransToReal(C2star));
  D2C2values.push_back(TransToReal(D2shepp));
  D2C2values.push_back(TransToReal(C2shepp));
  //cout << D2C2values[0] << endl;

  return D2C2values;
}

string D2StatNames[6] = { "D2NGS", "d2NGS", "D2starNGS", "d2starNGS", "D2sheppNGS", "d2sheppNGS" };


///////////////////////////////////////////////////////////////////////////////////////////
//KmerCount.out [k] [sample data file] 
int main(int argc, char **argv)   //EDIT main(int argc, char *argv[])
{
	const double Sim_Thres = 0.9;
	int k = 4;
	string kStr = "4";
	int order = 0;

	if (argc >= 7)
	{
		std::string kmerDir(argv[1]);
		std::string contigsNamelistURL(argv[2]);
		std::string d2starOutputURL(argv[3]);
		std::string d2sheppOutputURL(argv[4]);
		std::string jobIDStr(argv[5]); int jobID = atoi(jobIDStr.c_str());
		std::string totalJobNumStr(argv[6]); int totalJobNum = atoi(totalJobNumStr.c_str());

		ofstream d2starOfs(d2starOutputURL.c_str(), std::ofstream::out);
		ofstream d2sheppOfs(d2sheppOutputURL.c_str(), std::ofstream::out);

		string line;
		vector<string> contigsNameList;
		ifstream contigStream(contigsNamelistURL.c_str());
		if (contigStream.is_open())
		{
			while (getline(contigStream, line))
			{
				if (!line.empty() && line != "")
					contigsNameList.push_back(line);
			}
			contigStream.close();
		}

		std::cout << contigsNameList.size() << std::endl;
		for (int i = 0; i < contigsNameList.size(); ++i)
		{
			if ((i % totalJobNum) != jobID)
			{
				std::cout << "skip the contig id: " << i << std::endl;
				continue;
			}
			KmerInfo* speciesKmerInfoA = new KmerInfo();

			std::string kmerFilePathNameA = kmerDir + contigsNameList[i] + "_k" + kStr + "_singleStrand_wordcount";
			loadKmerCountHash(kmerFilePathNameA, speciesKmerInfoA, "kmerCount");
			loadKmerCountHash(kmerFilePathNameA, speciesKmerInfoA, "kmerOrder+1");
			pwIID(ZI, k, speciesKmerInfoA);

			for (int j = i + 1; j < contigsNameList.size(); ++j)
			{
				KmerInfo* speciesKmerInfoB = new KmerInfo();

				std::string kmerFilePathNameB = kmerDir + contigsNameList[j] + "_k" + kStr + "_singleStrand_wordcount";
				loadKmerCountHash(kmerFilePathNameB, speciesKmerInfoB, "kmerCount");
				loadKmerCountHash(kmerFilePathNameB, speciesKmerInfoB, "kmerOrder+1");
				pwIID(ZI, k, speciesKmerInfoB);

				vector<double> D2C2NGSvalues = D2C2computeNGS(ZI, k, speciesKmerInfoA, speciesKmerInfoB);

				/*std::cout << "compute " << contigsNameList[i] << " and " << contigsNameList[j] << std::endl;
				for (int m = 0; m < D2C2NGSvalues.size(); m++)
				{
					std::cout << D2StatNames[m] << ", " << D2C2NGSvalues[m] << std::endl;
				}*/

				if (D2C2NGSvalues[3] >= Sim_Thres)
					d2starOfs << (i + 1) << "\t" << (j + 1) << "\t" << D2C2NGSvalues[3] << endl;

				if (D2C2NGSvalues[5] >= Sim_Thres)
					d2sheppOfs << (i + 1) << "\t" << (j + 1) << "\t" << D2C2NGSvalues[5] << endl;

				delete speciesKmerInfoB;
			}
			std::cout << "Finished " << contigsNameList[i] << std::endl;

			delete speciesKmerInfoA;
		}

		d2starOfs.close();
		d2sheppOfs.close();
	}

  return 0;
}

