
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

//#include<ext/hash_map>
//using __gnu_cxx::hash_map;

using namespace std;



struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};


int ZI = 4;
//int k = 3;
//int order = 1;


//receptacle of kmer counts
unordered_map<unsigned long,unsigned long> HashTable[2];
// kmer length = order

unordered_map<unsigned long,unsigned long> HashTableK_1[2];
// kmer length = k - 2
unordered_map<unsigned long,unsigned long> HashTableK_2[2];

unordered_map<unsigned long,unsigned long> HashTableOrder[2];
// kmer length = order + 1
unordered_map<unsigned long,unsigned long> HashTableOrder_1[2];

//receptacle of Pw(probability of a kmer word) 
unordered_map<unsigned long,SCIENTIFIC_NUMBER > HashPw[2];

unsigned long totalKmer[2];
unsigned long totalK_1[2];
unsigned long totalK_2[2];
unsigned long totalOrder[2];
unsigned long totalOrder_1[2];





// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
  double dReal=0;
  dReal = dSci.value * pow(10,dSci.factor);
  return dReal;
  
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


void loadKmerCountHash(char *currentKmerFilePathName, int speciesID, string kmerUsage)
{
  ifstream fin(currentKmerFilePathName);  
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
      HashTable[speciesID][currentKmerID] = currentKmerCount;
      totalKmer[speciesID] = totalKmer[speciesID] + currentKmerCount;
      //cout << "kmerCount " << currentKmerCount << endl;
    }else if(kmerUsage == "kmerCount_1"){
      HashTableK_1[speciesID][currentKmerID] = currentKmerCount;
      totalK_1[speciesID] = totalK_1[speciesID] + currentKmerCount;
      //cout << "kmerOrder " << currentKmerCount << endl;
    }else if(kmerUsage == "kmerCount_2"){
      HashTableK_2[speciesID][currentKmerID] = currentKmerCount;
      totalK_2[speciesID] = totalK_2[speciesID] + currentKmerCount;
      //cout << "kmerOrder " << currentKmerCount << endl;
    }else if(kmerUsage == "kmerOrder"){
      HashTableOrder[speciesID][currentKmerID] = currentKmerCount;
      totalOrder[speciesID] = totalOrder[speciesID] + currentKmerCount;
      //cout << "kmerOrder " << currentKmerCount << endl;
    }else if(kmerUsage == "kmerOrder+1"){
      HashTableOrder_1[speciesID][currentKmerID] = currentKmerCount;
      totalOrder_1[speciesID] = totalOrder_1[speciesID] + currentKmerCount;
      //cout << "kmerOrder+1 " << currentKmerCount << endl;
    }else{cout << "ERROR: wrong kmerUsage" << endl;}
    
  }
  
  //for(int key=0; key<pow(ZI,k); key++ )
  //{
  //cout << key << " " << HashTable[speciesID][key] << endl;
  //}
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
  


void pwIID(int ZI, int k, int speciesID)
{
  // compute pw for each kmer word
  unordered_map<unsigned long,SCIENTIFIC_NUMBER > probIID;
  //cout << HashTableOrder_1[speciesID][0]/double(totalOrder_1[speciesID]) << endl;
  for(int index = 0; index < ZI; index++)
  { 
    //cout << HashTableOrder_1[speciesID][index] << endl;
    probIID[index] = TransToScientific(HashTableOrder_1[speciesID][index]/double(totalOrder_1[speciesID]));
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
    
    HashPw[speciesID][ten] = pw;
  }
}



void pwMC(int ZI, int k, int order, int speciesID, vector<SCIENTIFIC_NUMBER>& iniProb, vector< vector<SCIENTIFIC_NUMBER> > &transMatrix)
{
  // compute the transition probability matrix
  int transColSize = ZI;
  int transRowSize = pow(ZI, order);
  //vector< vector<SCIENTIFIC_NUMBER> > transMatrix(transRowSize, vector<SCIENTIFIC_NUMBER>(transColSize));
  //vector<SCIENTIFIC_NUMBER> iniProb;
  for(int currentRow = 0; currentRow < transRowSize; currentRow++)
  {
    vector<int> currentRowFour = ten2four(currentRow, order);
    unsigned long countBelow = HashTableOrder[speciesID][currentRow];
    double probBelow = countBelow/double(totalOrder[speciesID]);
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
      unsigned long countAbove = HashTableOrder_1[speciesID][currentRowColTen];
      double probAbove = countAbove/double(totalOrder_1[speciesID]);
      SCIENTIFIC_NUMBER probAboveSci = TransToScientific(probAbove);

      if( probBelow != 0 )
      {
        transMatrix[currentRow][currentCol] = SciMultiple(probAboveSci,inv_probBelowSci);
      }else{
        
        transMatrix[currentRow][currentCol] = TransToScientific(0);
      }
      
      rowSum = SciAddition(rowSum, transMatrix[currentRow][currentCol]);
    
    //cout << currentRow << ", " << currentCol << ", " << countAbove << ", " << countBelow << ", " << totalOrder_1[speciesID] << ", " << totalOrder[speciesID] << ", " << countAbove/double(countBelow) << ", " << probAbove << ", " << probBelow << ", " << TransToReal(transMatrix[currentRow][currentCol]) << endl;
    //cout << totalKmer[speciesID] << ", " << totalOrder[speciesID] << ", "<< totalOrder_1[speciesID] << endl;
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
    
    HashPw[speciesID][currentKmerTen] = pw;
    //printFour(currentKmerFour);
    //cout << " pw:" << TransToReal(pw) << endl;
  }
  
}



vector<double> D2C2computeLongseq(int ZI, int k)
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
    //vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
    
    //unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
    //cout << currentKmerRevTen << endl;
    //printFour(currentKmerFour);
    //printFour(currentKmerFour);
    //printFour(currentKmerRevFour);
    
    SCIENTIFIC_NUMBER p_w[2], EX_w[2], X_w[2], X_w_tilde[2], X_w_tilde_var[2];
    
    for(int speciesID=0; speciesID<2; speciesID++)
    { 
      p_w[speciesID] = HashPw[speciesID][currentKmerTen];
      EX_w[speciesID] = SciMultiple(p_w[speciesID], totalKmer[speciesID]);
      X_w[speciesID] = TransToScientific(HashTable[speciesID][currentKmerTen]);
      //cout << "Species" << speciesID << " Ten:" << currentKmerTen << " Four: " ;
      //printFour(currentKmerFour);
      //cout << "forward:" << HashTable[speciesID][currentKmerTen] << " reverse:" << HashTable[speciesID][currentKmerRevTen] << TransToReal(X_w[speciesID]) << endl;
      X_w_tilde[speciesID] = SciAddition(X_w[speciesID], SciNegative(EX_w[speciesID]));
      
      //cout << TransToReal(p_w[speciesID]) << ", " << TransToReal(EX_w[speciesID]) << ", " << TransToReal(X_w[speciesID]) << endl;
      
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






vector<double> D2C2computeNGS(int ZI, int k)
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
      p_w[speciesID] = SciAddition(HashPw[speciesID][currentKmerTen], HashPw[speciesID][currentKmerRevTen]);
      EX_w[speciesID] = SciMultiple(p_w[speciesID], totalKmer[speciesID]);
      X_w[speciesID] = TransToScientific(HashTable[speciesID][currentKmerTen] + HashTable[speciesID][currentKmerRevTen]);
      //cout << "Species" << speciesID << " Ten:" << currentKmerTen << " Four: " ;
      //printFour(currentKmerFour);
      //cout << "forward:" << HashTable[speciesID][currentKmerTen] << " reverse:" << HashTable[speciesID][currentKmerRevTen] << TransToReal(X_w[speciesID]) << endl;
      X_w_tilde[speciesID] = SciAddition(X_w[speciesID], SciNegative(EX_w[speciesID]));
      
      //cout << TransToReal(p_w[speciesID]) << ", " << TransToReal(EX_w[speciesID]) << ", " << TransToReal(X_w[speciesID]) << endl;
      
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




vector<double> EuMaChDistNGS(int ZI, int k)
{
  //cout << TransToReal(D2) << endl;
  SCIENTIFIC_NUMBER EuDist ; EuDist.value=0; EuDist.factor=0;
  SCIENTIFIC_NUMBER MaDist ; MaDist.value=0; MaDist.factor=0;
  SCIENTIFIC_NUMBER ChDist ; ChDist.value=0; ChDist.factor=0;
  
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
    unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
    
    // frequency
    SCIENTIFIC_NUMBER X1 = TransToScientific( ( HashTable[0][currentKmerTen] + HashTable[0][currentKmerRevTen] ) / double ( 2 * totalKmer[0] ) );
    SCIENTIFIC_NUMBER X2 = TransToScientific( ( HashTable[1][currentKmerTen] + HashTable[1][currentKmerRevTen] ) / double ( 2 * totalKmer[1] ) );
    
    SCIENTIFIC_NUMBER diff = SciAddition(X1, SciNegative(X2)) ;
    diff.value = fabs(diff.value);
    
    EuDist = SciAddition( EuDist, SciPow( diff , 2 ) );
    
    MaDist = SciAddition( MaDist, diff ) ;
    
    if( TransToReal(diff) > TransToReal(ChDist) )
    {
      ChDist = diff;
    }
  }
  
  vector<double> EuMaChDistvalues;
  EuMaChDistvalues.push_back(TransToReal(SciPow(EuDist, 0.5)));
  EuMaChDistvalues.push_back(TransToReal(MaDist));
  EuMaChDistvalues.push_back(TransToReal(ChDist));

  return EuMaChDistvalues;
  
}




vector<double> modifiedEuDistNGS(int ZI, int k)
{
	//cout << TransToReal(D2) << endl;
	// expected under MC model, divide
	SCIENTIFIC_NUMBER EuFDist ; EuFDist.value=0; EuFDist.factor=0;
	
	for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
	{
		vector<int> currentKmerFour = ten2four(currentKmerTen, k);
		vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
		unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
		
		// frequency
		SCIENTIFIC_NUMBER Nw1 = TransToScientific( HashTable[0][currentKmerTen] + HashTable[0][currentKmerRevTen] ) ;
		SCIENTIFIC_NUMBER pw1 = SciAddition( HashPw[0][currentKmerTen], HashPw[0][currentKmerRevTen]);
		SCIENTIFIC_NUMBER Ew1 = SciMultiple( pw1, totalKmer[0]);
		SCIENTIFIC_NUMBER Fw1;
		Fw1.value = 0; Fw1.factor = 0;
		if (Ew1.value != 0)
		{
			Fw1 = SciMultiple( Nw1, SciInverse(Ew1) );
		}
		
		SCIENTIFIC_NUMBER Nw2 = TransToScientific( HashTable[1][currentKmerTen] + HashTable[1][currentKmerRevTen] ) ;
		SCIENTIFIC_NUMBER pw2 = SciAddition( HashPw[1][currentKmerTen], HashPw[1][currentKmerRevTen] );
		SCIENTIFIC_NUMBER Ew2 = SciMultiple( pw2, totalKmer[1]);
		SCIENTIFIC_NUMBER Fw2;
		Fw2.value = 0; Fw2.factor = 0;
		if (Ew2.value != 0)
		{
			Fw2 = SciMultiple( Nw2, SciInverse(Ew2));
		}
		
		//cout << "Species" << speciesID << " Ten:" << currentKmerTen << " Four: " ;
		//printFour(currentKmerFour);
		//cout << "forward:" << HashTable[speciesID][currentKmerTen] << " reverse:" << HashTable[speciesID][currentKmerRevTen] << TransToReal(X_w[speciesID]) << endl;
		
		SCIENTIFIC_NUMBER diffF = SciAddition(Fw1, SciNegative(Fw2)) ;
		
		EuFDist = SciAddition( EuFDist, SciPow( diffF , 2 ) );

		
		//cout << TransToReal(Nw1) << ", " << TransToReal(Ew1) << ", " << TransToReal(Fw1) << ", " << TransToReal(Cw1) << ", "  << TransToReal(diffF) << endl;
		//cout << TransToReal(Nw2) << ", " << TransToReal(Ew2) << ", " << TransToReal(Fw2) << ", " << TransToReal(Cw2) << ", " << TransToReal(diffC) << endl;
		
	}
	
	vector<double> EuFDistvalues;
	EuFDistvalues.push_back(TransToReal(SciPow(EuFDist, 0.5)));

	return EuFDistvalues;
	
}






// Willner et al. Di, Tri, Tetra
double WillnerDiNGS(int ZI, int k)
{
  
  double deltaDi = 0;
  //cout << "test" << ", k" << k << endl;

  // Di
  if( k == 2 )
  {

    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    {
      //cout << "test" << endl;
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
  //    cout << "test" << endl;
      
      // frequency
      double fab1 = ( HashTable[0][currentKmerTen] + HashTable[0][currentKmerRevTen] ) / double( 2 * totalKmer[0] );
      double fab2 = ( HashTable[1][currentKmerTen] + HashTable[1][currentKmerRevTen] ) / double( 2 * totalKmer[1] );
      
      //cout << "fab1" << fab1 << endl;
      
      double fa1 = ( HashTableOrder_1[0][currentKmerFour[0]] + HashTableOrder_1[0][3-currentKmerFour[0]] ) / double( 2 * totalOrder_1[0] );
      double fb1 = ( HashTableOrder_1[0][currentKmerFour[1]] + HashTableOrder_1[0][3-currentKmerFour[1]] ) / double( 2 * totalOrder_1[0] );
      
      double fa2 = ( HashTableOrder_1[1][currentKmerFour[0]] + HashTableOrder_1[1][3-currentKmerFour[0]] ) / double( 2 * totalOrder_1[1] );
      double fb2 = ( HashTableOrder_1[1][currentKmerFour[1]] + HashTableOrder_1[1][3-currentKmerFour[1]] ) / double( 2 * totalOrder_1[1] );
      
      deltaDi = deltaDi + fabs( fab1 / ( fa1 * fb1) - fab2 / ( fa2 * fb2) );
      
    } 
    return deltaDi;
  }
  
  double gammaTri = 0;
  if( k == 3 )
  { 
    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    {
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
      // frequency
      double fabc1 = ( HashTable[0][currentKmerTen] + HashTable[0][currentKmerRevTen] ) / double( 2 * totalKmer[0] ) ;
      double fabc2 = ( HashTable[1][currentKmerTen] + HashTable[1][currentKmerRevTen] ) / double( 2 * totalKmer[1] ) ;
      
      vector<double> fTwoLetters1; // 3 in total
      vector<double> fTwoLetters2; // 3 in total
      for( int missPos = 0; missPos < 4; missPos++)
      {
        double wordFreq1=0;
        double wordFreq2=0;
        for( int missLetter = 0; missLetter < 4; missLetter++)
        {
          vector<int> currentWordFour = currentKmerFour;
          currentWordFour[missPos] = missLetter;
          unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
          unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
          wordFreq1 = wordFreq1 + ( HashTable[0][currentWordTen] + HashTable[0][currentWordRevTen] ) ;
          wordFreq2 = wordFreq2 + ( HashTable[1][currentWordTen] + HashTable[1][currentWordRevTen] ) ;
        }
        fTwoLetters1.push_back( wordFreq1 / double (2 * totalKmer[0] ) ) ;
        fTwoLetters2.push_back( wordFreq2 / double (2 * totalKmer[1] ) ) ;

      }
      
      vector<double> fOneLetters1; // 3 in total
      vector<double> fOneLetters2; // 3 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          double wordFreq1=0;
          double wordFreq2=0;
          for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
          {
            for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
            {
              //cout << "missPosFirst," << missPosFirst << ",missPosSecond," << missPosSecond << endl;
              //cout << "missLetterFirst," << missLetterFirst << ",missLetterSecond," << missLetterSecond << endl; 
              vector<int> currentWordFour = currentKmerFour;
              currentWordFour[missPosFirst] = missLetterFirst;
              currentWordFour[missPosSecond] = missLetterSecond;
              unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
              unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
              wordFreq1 = wordFreq1 + ( HashTable[0][currentWordTen] + HashTable[0][currentWordRevTen] ) ;
              wordFreq2 = wordFreq2 + ( HashTable[1][currentWordTen] + HashTable[1][currentWordRevTen] ) ;
              //cout << "missPosFirst," << missPosFirst << ",missPosSecond," << missPosSecond << endl; 
              //printFour(currentKmerFour);
              //cout <<  "wordFreq," << wordFreq1 << endl;
            }
          }
          fOneLetters1.push_back( wordFreq1 / double (2 * totalKmer[0] ) ) ;
          fOneLetters2.push_back( wordFreq2 / double (2 * totalKmer[1] ) ) ;
        }
      }
      //cout << fOneLetters1[0] << "," << fOneLetters1[1] << "," << fOneLetters1[2] << endl; 
      //cout << fOneLetters1[0] << "," << fOneLetters1[1] << "," << fOneLetters1[2] << endl;
      double gamma1 = fabc1 * fOneLetters1[0] * fOneLetters1[1] * fOneLetters1[2] / ( fTwoLetters1[0] * fTwoLetters1[1] * fTwoLetters1[2] ) ;
      double gamma2 = fabc2 * fOneLetters2[0] * fOneLetters2[1] * fOneLetters2[2] / ( fTwoLetters2[0] * fTwoLetters2[1] * fTwoLetters2[2] ) ;
      gammaTri = gammaTri + fabs( gamma1 - gamma2) ;
      //cout << gamma1 << ","  << gamma2 << endl;
      //cout << fabs(gamma1-gamma2) << ".." << gammaTri << endl;
    }
    //cout << a << "," << gammaTri << endl;
    return gammaTri;
  }
  
  double tauTetra = 0;
  if( k == 4)
  {
    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    { 
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
      double fabcd1 = ( HashTable[0][currentKmerTen] + HashTable[0][currentKmerRevTen] ) / double( 2 * totalKmer[0] ) ;
      double fabcd2 = ( HashTable[1][currentKmerTen] + HashTable[1][currentKmerRevTen] ) / double( 2 * totalKmer[1] ) ;
      
      vector<double> fThreeLetters1; // 4 in total
      vector<double> fThreeLetters2; // 4 in total
      for( int missPos = 0; missPos < 4; missPos++)
      {
        double wordFreq1=0;
        double wordFreq2=0;
        for( int missLetter = 0; missLetter < 4; missLetter++)
        {
          vector<int> currentWordFour = currentKmerFour;
          currentWordFour[missPos] = missLetter;
          unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
          unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
          wordFreq1 = wordFreq1 + ( HashTable[0][currentWordTen] + HashTable[0][currentWordRevTen] ) ;
          wordFreq2 = wordFreq2 + ( HashTable[1][currentWordTen] + HashTable[1][currentWordRevTen] ) ;
        }
        fThreeLetters1.push_back( wordFreq1 / double (2 * totalKmer[0] ) ) ;
        fThreeLetters2.push_back( wordFreq2 / double (2 * totalKmer[1] ) ) ;
        
      }
      
      vector<double> fTwoLetters1; // 6 in total
      vector<double> fTwoLetters2; // 6 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          double wordFreq1=0;
          double wordFreq2=0;
          for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
          {
            for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
            {
              vector<int> currentWordFour = currentKmerFour;
              currentWordFour[missPosFirst] = missLetterFirst;
              currentWordFour[missPosSecond] = missLetterSecond;
              unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
              unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
              wordFreq1 = wordFreq1 + ( HashTable[0][currentWordTen] + HashTable[0][currentWordRevTen] ) ;
              wordFreq2 = wordFreq2 + ( HashTable[1][currentWordTen] + HashTable[1][currentWordRevTen] ) ;
            }
          }
          fTwoLetters1.push_back( wordFreq1 / double (2 * totalKmer[0] ) ) ;
          fTwoLetters2.push_back( wordFreq2 / double (2 * totalKmer[1] ) ) ;
        }
      }

      vector<double> fOneLetters1; // 4 in total
      vector<double> fOneLetters2; // 4 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          for( int missPosThird = missPosSecond + 1; missPosThird < 4; missPosThird++)
          {
            double wordFreq1=0;
            double wordFreq2=0;
            for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
            {
              for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
              {
                for( int missLetterThird = 0; missLetterThird < 4; missLetterThird++)
                {
                  vector<int> currentWordFour = currentKmerFour;
                  currentWordFour[missPosFirst] = missLetterFirst;
                  currentWordFour[missPosSecond] = missLetterSecond;
                  currentWordFour[missPosThird] = missLetterThird;
                  unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
                  unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
                  wordFreq1 = wordFreq1 + ( HashTable[0][currentWordTen] + HashTable[0][currentWordRevTen] ) ;
                  wordFreq2 = wordFreq2 + ( HashTable[1][currentWordTen] + HashTable[1][currentWordRevTen] ) ;
                }
              }
            }
            
            fOneLetters1.push_back( wordFreq1 / double (2 * totalKmer[0] ) ) ;
            fOneLetters2.push_back( wordFreq2 / double (2 * totalKmer[1] ) ) ;
            
            //cout << missPosFirst << "," << missPosSecond << "," << missPosThird << endl;
            //cout << fOneLetters1.size() << "," << fOneLetters1[fOneLetters1.size()-1] << endl;
          }
        }
      }
      //cout << "fTwoLetters1_0," << fTwoLetters1[0] << endl;
      //cout << fOneLetters1[1] << endl;
      //cout << fTwoLetters1[1] << endl;
      //cout << fThreeLetters1[1] << endl;
      double tau1 = fabcd1 * fTwoLetters1[0] * fTwoLetters1[1] * fTwoLetters1[2] * fTwoLetters1[3] * fTwoLetters1[4] * fTwoLetters1[5] / ( fThreeLetters1[0] * fThreeLetters1[1] * fThreeLetters1[2] * fThreeLetters1[3] * fOneLetters1[0] * fOneLetters1[1] * fOneLetters1[2] * fOneLetters1[3]) ;
      
      double tau2 = fabcd2 * fTwoLetters2[0] * fTwoLetters2[1] * fTwoLetters2[2] * fTwoLetters2[3] * fTwoLetters2[4] * fTwoLetters2[5] / ( fThreeLetters2[0] * fThreeLetters2[1] * fThreeLetters2[2] * fThreeLetters2[3] * fOneLetters2[0] * fOneLetters2[1] * fOneLetters2[2] * fOneLetters2[3]) ;
      
      //cout << "tau1," << tau1 << endl;
      
      tauTetra = tauTetra + fabs(tau1 - tau2) ;
    }

    return tauTetra ;
  }
  
  
}




/// need to write it as complimentary chains
double HAOcompute(int ZI, int k)
{
  SCIENTIFIC_NUMBER HAO_above;
  HAO_above.value = 0; HAO_above.factor = 0;
  SCIENTIFIC_NUMBER HAO_below[2]; 
  HAO_below[0].value = 0; HAO_below[0].factor = 0; 
  HAO_below[1].value = 0; HAO_below[1].factor = 0; 
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++ )
  {
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << endl << endl;
    //printFour(currentKmerFour);
    vector<int> currentKmerRemoveLastFour;
    vector<int> currentKmerRemoveFirstFour;
    vector<int> currentKmerRemoveTwoFour;
    for(int position = 0; position < k; position++)
    {
      if( position != k)
      {
        currentKmerRemoveLastFour.push_back(currentKmerFour[position]);
      }
      if( position != 0)
      {
        currentKmerRemoveFirstFour.push_back(currentKmerFour[position]);
      }
      if( position != 0 && position != k)
      {
        currentKmerRemoveTwoFour.push_back(currentKmerFour[position]);
      }
    }
    unsigned long currentKmerRemoveLastTen = four2ten(currentKmerRemoveLastFour, (k-1));
    unsigned long currentKmerRemoveFirstTen = four2ten(currentKmerRemoveFirstFour, (k-1));
    unsigned long currentKmerRemoveTwoTen = four2ten(currentKmerRemoveTwoFour, (k-2));
    
    SCIENTIFIC_NUMBER avalue[2]; 
    avalue[0].value = 0; avalue[0].factor = 0; 
    avalue[1].value = 0; avalue[1].factor = 0;
    for(int speciesID=0; speciesID<2; speciesID++)
    { 
      SCIENTIFIC_NUMBER pw0; pw0.value = 0; pw0.factor = 0;
      if(HashTableK_2[speciesID][currentKmerRemoveTwoTen] != 0)
      {
        //cout << "HashTableK_1 " << "species " << speciesID << " word " << currentKmerRemoveLastTen << " kmercount " << HashTableK_1[speciesID][currentKmerRemoveLastTen] << " totalK_1 " << totalK_1[speciesID] << "; " << HashTableK_1[speciesID][currentKmerRemoveFirstTen] << " totalK_1 " << totalK_1[speciesID] << "; " << HashTableK_2[speciesID][currentKmerRemoveTwoTen] << " totalK_2 " << totalK_2[speciesID] << endl;
        SCIENTIFIC_NUMBER pw0_up1 = TransToScientific(HashTableK_1[speciesID][currentKmerRemoveLastTen]/double(totalK_1[speciesID]));
        SCIENTIFIC_NUMBER pw0_up2 = TransToScientific(HashTableK_1[speciesID][currentKmerRemoveFirstTen]/double(totalK_1[speciesID]));
        SCIENTIFIC_NUMBER pw0_down = TransToScientific(HashTableK_2[speciesID][currentKmerRemoveTwoTen]/double(totalK_2[speciesID]));
        
        pw0 = SciMultiple(SciMultiple(pw0_up1, pw0_up2),SciInverse(pw0_down));
        //cout << HashTableK_1[speciesID][currentKmerRemoveLastTen]/totalK_1[speciesID] << " " << TransToReal(pw0_up1) << " " << TransToReal(pw0_up2) << " " << TransToReal(pw0_down) << " " << TransToReal(pw0) << endl;
        //fw0 = totalKmer[speciesID] * pw0;
        
      }
      
      SCIENTIFIC_NUMBER pw = TransToScientific(HashTable[speciesID][currentKmerTen]/double(totalKmer[speciesID]));
      //cout << "pw " << TransToReal(pw) << " pw0 " << TransToReal(pw0) << endl;
      if(pw0.value != 0)
      {
        avalue[speciesID] = SciMultiple(SciAddition(pw, SciNegative(pw0)), SciInverse(pw0));
      }
      //cout << TransToReal(pw) << " " << TransToReal(pw0) << endl;
      //cout << TransToReal(avalue[speciesID]) << endl;
      HAO_below[speciesID] = SciAddition(HAO_below[speciesID], SciPow(avalue[speciesID],2));
    }
    
    HAO_above = SciAddition(HAO_above, SciMultiple(avalue[0], avalue[1]));
    
  }
  //cout << TransToReal(HAO_above) << " " << TransToReal(HAO_below[0]) << " " << TransToReal(HAO_below[1]) << endl;
  
  SCIENTIFIC_NUMBER HAOvalueSCI = SciMultiple(HAO_above, SciInverse(SciPow(SciMultiple(HAO_below[0], HAO_below[1]) , 0.5)));                                                                    
  double HAOvalue = TransToReal(HAOvalueSCI);     
  //cout << HAOvalue << endl;
  
  //HAOvalue = (1 - HAOvalue)/2;
  //cout << HAOvalue << endl;
  return HAOvalue;
}






double JScomputeIID(vector<SCIENTIFIC_NUMBER>& iniProb)
{
  int size = iniProb.size();
  SCIENTIFIC_NUMBER entropy;
  entropy.value = 0; entropy.factor = 0;
  for(int currentIndex = 0; currentIndex < size; currentIndex++)
  {
    // compute JS distance
    if( TransToReal(iniProb[currentIndex]) != 0 ){
      
      SCIENTIFIC_NUMBER logTransProbSci = TransToScientific(log2(TransToReal(iniProb[currentIndex])));
      entropy = SciAddition(entropy, SciMultiple(iniProb[currentIndex], logTransProbSci));
    }else{
      //cout << "trans prob = 0, no addition on the entropy." << endl;
    }
  }
  return -TransToReal(entropy);
}






double JScomputeMC(vector<SCIENTIFIC_NUMBER>& iniProb, vector< vector<SCIENTIFIC_NUMBER> >& transMatrix)
{
  int transRowSize = iniProb.size();
  int transColSize = transMatrix[0].size();
  SCIENTIFIC_NUMBER entropyOverRow;
  entropyOverRow.value = 0; entropyOverRow.factor = 0;
  for(int currentRow = 0; currentRow < transRowSize; currentRow++)
  {
    // normalize trans matrix
    SCIENTIFIC_NUMBER entropyOverCol;
    entropyOverCol.value = 0; entropyOverCol.factor = 0;
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      // compute JS distance
      if( TransToReal(transMatrix[currentRow][currentCol]) != 0 ){
        
        SCIENTIFIC_NUMBER logTransProbSci = TransToScientific(log2(TransToReal(transMatrix[currentRow][currentCol])));
        entropyOverCol = SciAddition(entropyOverCol, SciMultiple(transMatrix[currentRow][currentCol], logTransProbSci));
      }else{
        //cout << "trans prob = 0, no addition on the entropy." << endl;
      }
    }
    
    entropyOverRow = SciAddition(entropyOverRow, SciMultiple(iniProb[currentRow], entropyOverCol));
    
  }
  return -TransToReal(entropyOverRow);
  
}






// not necessary to make it to include complimentary
double S2compute(int ZI, int k)
{

  //cout << TransToReal(D2) << endl;
  SCIENTIFIC_NUMBER S2;
  S2.factor = 0; S2.value = 0;
  
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << currentKmerRevTen << endl;
    //printFour(currentKmerFour);
    
    SCIENTIFIC_NUMBER phi[2];
    
    for(int speciesID=0; speciesID<2; speciesID++)
    { 
      phi[speciesID] = SciMultiple(TransToScientific( HashTable[speciesID][currentKmerTen] / double(totalKmer[speciesID]) ), HashPw[speciesID][currentKmerTen] ) ;
      //cout << HashTable[speciesID][currentKmerTen] << "," << totalKmer[speciesID] << ", pw: " << TransToReal(HashPw[speciesID][currentKmerTen]) << endl;
      
      //cout << pow(ZI, k) << ", " << currentKmerTen << ", phi: " << TransToReal(phi[speciesID]) << endl;

    }

    
    if( phi[0].value != 0 )
    {
      SCIENTIFIC_NUMBER tmp1 = SciLn(SciMultiple(phi[0], 2)) ;
      
      SCIENTIFIC_NUMBER tmp2 = SciLn(SciAddition(phi[0], phi[1])) ;
      
      SCIENTIFIC_NUMBER tmp3 = SciMultiple(phi[0], SciAddition( tmp1, SciNegative(tmp2))) ;
      tmp3.value = fabs(tmp3.value) ;
      
      //cout << "tmp1: " << TransToReal(tmp1) << ", tmp2: " << TransToReal(tmp2) << ", tmp1-tmp2: " << TransToReal(SciAddition( tmp1, SciNegative(tmp2))) << ", tmp3: " << TransToReal(tmp3) << endl;
      
      S2 = SciAddition( S2, tmp3 ) ;
      
    }
    
    if( phi[1].value != 0 )
    {
      
      SCIENTIFIC_NUMBER tmp4 = SciLn(SciMultiple(phi[1], 2)) ;
      SCIENTIFIC_NUMBER tmp2 = SciLn(SciAddition(phi[0], phi[1])) ;
      SCIENTIFIC_NUMBER tmp5 = SciMultiple(phi[1], SciAddition( tmp4, SciNegative(tmp2) )) ;
      tmp5.value = fabs(tmp5.value) ;
      
      S2 = SciAddition( S2, tmp5 ) ;
      
    }
    
    //cout << pow(ZI, k) << ", " << currentKmerTen << ", " << TransToReal(phi[0]) << ", " << TransToReal(phi[1]) << ", " << TransToReal(S2) << endl;
  }
  
  // so far S2 is very small, after dividing by pow(ZI, k), this term is almost 0!!!
  // so S2 is almost always 2*ln2!!! I dont get it!
  
  S2 = SciMultiple( S2, SciInverse(TransToScientific(pow(ZI, k))) ) ;
  S2 = SciAddition( S2, TransToScientific(2*log(2)) ) ;
  
  return TransToReal(S2);

}








                            


///////////////////////////////////////////////////////////////////////////////////////////
//KmerCount.out [k] [sample data file] 
int main(int argc, char **argv)   //EDIT main(int argc, char *argv[])
{
	//bool doubleStrand = false;
	char *speciesName[2];
	vector<int> order(2,0);
	
	int k = 0;
	char kstr[5];
	
	char *inputFileDir[2];
	inputFileDir[0] = NULL;
	inputFileDir[1] = NULL;
	char *outputFileDir = NULL;
	
	
	/* getOptions from command line */
	int c;
	
	while (1)
	{
		static struct option long_options[] =
		{
			//{"doubleStrand",     no_argument,       0, 'd'},
			// necessary arguments: kvalue and inputFileName !!
			{"species1",  required_argument, 0, 'a'},
			{"order1",  required_argument, 0, 'b'},
			{"species2",  required_argument, 0, 'c'},
			{"order2",  required_argument, 0, 'd'},
			{"kvalue",  required_argument, 0, 'k'},
			{"inputFileDir1",  required_argument, 0, 'i'},
			{"inputFileDir2",  required_argument, 0, 'j'},
			{"outputFileDir",  required_argument, 0, 'o'},
			{0,         0,                 0,  0 }
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		c = getopt_long (argc, argv, "a:b:c:d:k:i:j:o:",
										 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c)
		{

			//case 'd':
			//	doubleStrand = true;
			//	puts ("option -d: consider complimentary word. ");
			//	break;
				
			/* input kmer files should be singleStrand count! */
			 
			case 'a':
				speciesName[0] = optarg;
				//printf ("option -a, name of 1st species, with value `%s'\n", optarg);
				break;
				
			case 'b':
				order[0] = atoi(optarg);
				//printf ("option -b, order of 1st species, with value `%s'\n", optarg);
				break;
				
			case 'c':
				speciesName[1] = optarg;
				//printf ("option -c, name of 2nd species, with value `%s'\n", optarg);
				break;
				
			case 'd':
				order[1] = atoi(optarg);
				//printf ("option -d, order of 2nd species, with value `%s'\n", optarg);
				break;
				
			case 'k':
				k = atoi(optarg) ;
				sprintf(kstr, "%d", k);
				//printf ("option -k, the value of k, with value `%s'\n", optarg);
				break;
				
			case 'i':
				inputFileDir[0] = optarg;
				//printf ("option -i, the input kmer count file directory, with value `%s'\n", optarg);
				break;
				
			case 'j':
				inputFileDir[1] = optarg;
				//printf ("option -j, the input kmer count file directory, with value `%s'\n", optarg);
				break;
				
			case 'o':
				outputFileDir = optarg;
				//printf ("option -o, the output directory, with value `%s'\n", optarg);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				abort ();
		}
	}

	for(int speciesID=0; speciesID < 2; speciesID++)
	{
		if (inputFileDir[speciesID] != NULL) {
			std::string inputPath = inputFileDir[speciesID];
			if (inputPath[inputPath.length()-1] != '/') {
				strcat(inputFileDir[speciesID],"/");
			}
			//cout << "inputFileDir: " << inputFileDir[speciesID] << endl;
		}
	}
	
	if (outputFileDir != NULL) {
		std::string outputPath = outputFileDir;
		if (outputPath[outputPath.length()-1] != '/') {
			strcat(outputFileDir,"/");
		}
		//cout << "outputFileDir: " << outputFileDir << endl;
		
		// 0520 if outputFileDir not exist, mkdir
		struct stat st;
		if(stat(outputFileDir, &st) != 0)
		{
			//printf(" mkdir outputFileDir \n");
			char mkcmd[100];
			sprintf(mkcmd, "mkdir -p %s", outputFileDir);
			system(mkcmd);
		}
		
		//int status = mkdir (outputFileDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		//if( status == -1 )
		//{
		//	cout << "WARNING: outputFileDir already exists. " << endl;
		//}
		
	}
	
	
	
	///////////////////////////////////////////////////////////////////////
	//////////////////// load the kmerFiles and compute pwMC //////////////
	///////////////////////////////////////////////////////////////////////

  vector<SCIENTIFIC_NUMBER> iniProb[2];
  vector<vector<vector<SCIENTIFIC_NUMBER> > > transMatrix (2, vector<vector<SCIENTIFIC_NUMBER> > (pow(ZI,order[1]), vector<SCIENTIFIC_NUMBER> (ZI)));
  
  // load the kmer count hashtables of the two dataset

  for(int speciesID=0; speciesID<2; speciesID++)
  {
    //cout << "======= species " << speciesName[speciesID] << " order " << order[speciesID] << endl;
    char orderstr[5];
    sprintf(orderstr, "%d", order[speciesID]);
    char order1str[5];
    sprintf(order1str, "%d", (order[speciesID]+1));
    
    ///panfs/cmb-panasas2/renj/10species/rn5/rn5_d1_k14_d0_singleStrand_wordcount
    // load the kmer count
    //cout << "== load the kmer count ==" << endl;
    char kmerFilePathName[1000];
    strcpy(kmerFilePathName, inputFileDir[speciesID]);
    //strcat(kmerFilePathName, "kmerCount/");
    strcat(kmerFilePathName, speciesName[speciesID]); 
    strcat(kmerFilePathName, "_k");
    strcat(kmerFilePathName, kstr);
    strcat(kmerFilePathName, "_singleStrand_wordcount");
    //cout << "kmerFile: " << kmerFilePathName << endl;
    loadKmerCountHash(kmerFilePathName, speciesID, "kmerCount");
    
    if(order[speciesID] == 0)
    {
      // load the 1-kmer count for MC: k=order+1
      //cout << "== load the kmer count for MC, order is " << order[speciesID] << " ==" << endl;
      char kmerFileOrder1PathName[1000];
      strcpy(kmerFileOrder1PathName, inputFileDir[speciesID]);
      //strcat(kmerFileOrder1PathName, "/kmerCount/");
      strcat(kmerFileOrder1PathName, speciesName[speciesID]);
      strcat(kmerFileOrder1PathName, "_k");
      strcat(kmerFileOrder1PathName, order1str);
      strcat(kmerFileOrder1PathName, "_singleStrand_wordcount");
      //cout << "kmerFile: " << kmerFileOrder1PathName << endl;
      loadKmerCountHash(kmerFileOrder1PathName, speciesID, "kmerOrder+1");

      
      // compute the pw
      //cout << "== compute the pwIID for the words, for D2-type and S2 == " << endl;
      pwIID(ZI, k, speciesID);
        
    
    }else{
      // load the kmer count for MC: k=order
      //cout << "== load the kmer count for MC, order is " << orderstr << " ==" << endl;
      char kmerFileOrderPathName[1000];
      strcpy(kmerFileOrderPathName, inputFileDir[speciesID]);
      //strcat(kmerFileOrderPathName, "/kmerCount/");
      strcat(kmerFileOrderPathName, speciesName[speciesID]);
      strcat(kmerFileOrderPathName, "_k");
      strcat(kmerFileOrderPathName, orderstr);
      strcat(kmerFileOrderPathName, "_singleStrand_wordcount");
      //cout << "kmerFile: " << kmerFileOrderPathName << endl;
      loadKmerCountHash(kmerFileOrderPathName, speciesID, "kmerOrder");
      // load the kmer count for MC: k=order+1
      //cout << "== load the kmer count for MC, order+1 is " << order1str << " ==" << endl;
      char kmerFileOrder1PathName[1000];
      strcpy(kmerFileOrder1PathName, inputFileDir[speciesID]);
      //strcat(kmerFileOrder1PathName, "/kmerCount/");
      strcat(kmerFileOrder1PathName, speciesName[speciesID]);
      strcat(kmerFileOrder1PathName, "_k");
      strcat(kmerFileOrder1PathName, order1str);
      strcat(kmerFileOrder1PathName, "_singleStrand_wordcount");
      //cout << "kmerFile: " << kmerFileOrder1PathName << endl;
      loadKmerCountHash(kmerFileOrder1PathName, speciesID, "kmerOrder+1");

      // compute the pw
      //cout << "== compute the pw for the words == " << endl;
      pwMC(ZI, k, order[speciesID], speciesID, iniProb[speciesID], transMatrix[speciesID]);
      
      
    }
    
  }
  
  
  ///////////////////////////////////////////////////////////////////////  
  ////////////////////// computing the multiple statistics //////////////
  ///////////////////////////////////////////////////////////////////////
  char outputfile[1000];
  char orderstrA[5];
  sprintf(orderstrA, "%d", order[0]);
  char orderstrB[5];
  sprintf(orderstrB, "%d", order[1]);
  
  strcpy(outputfile, outputFileDir);
  strcat(outputfile, "k");
  strcat(outputfile,kstr);
  strcat(outputfile,"_");
  strcat(outputfile,speciesName[0]);
  strcat(outputfile,"_order");
  strcat(outputfile,orderstrA);
  strcat(outputfile,"_");
  strcat(outputfile,speciesName[1]);
  strcat(outputfile,"_order");
  strcat(outputfile,orderstrB);
  strcat(outputfile,"_statScores");
  //ofstream fout(outputfile);
  
	
  // 1. compute the D2 statistics
	string D2StatNames[12] = {"D2", "d2", "D2star", "d2star", "D2shepp", "d2shepp", "D2NGS", "d2NGS", "D2starNGS", "d2starNGS", "D2sheppNGS", "d2sheppNGS"} ;
  //cout << "== compute the D2 statistics, NOT consider complementary == " << endl;
  vector<double> D2C2Longseqvalues = D2C2computeLongseq(ZI, k);
  
  for(int i = 0; i < D2C2Longseqvalues.size(); i++)
  {
    //cout << i << endl;
    //fout << D2StatNames[i] << ", " << D2C2Longseqvalues[i] << endl;
    cout << D2StatNames[i] << ", " << D2C2Longseqvalues[i] << endl;
  }
  
  
  //cout << "== compute the D2 statistics, consider complementary == " << endl;
  vector<double> D2C2NGSvalues = D2C2computeNGS(ZI, k);
  //printFour(D2C2values);

  for(int i = 0; i < D2C2NGSvalues.size(); i++)
  {
    //cout << i << endl;
    //fout << D2StatNames[i+6] << ", " << D2C2NGSvalues[i] << endl ;
    cout << D2StatNames[i+6] << ", " << D2C2NGSvalues[i] << endl ;
  }

  
  
  // 2. compute the S2 statistic  
  //cout << "== compute the S2 statistics, considering single strand == " << endl;
  double S2 = S2compute(ZI, k);
  //fout << "S2, " << S2 << endl;
  cout << "S2, " << S2 << endl;
	
	
	
  // 3. compute HAO distance 
  //cout << "== compute the HAO distance, considering single strand == " << endl;
  if(k < 3){
    //cerr << "ERROR: There is no Hao distance for k < 3!" << endl;
    
  }else{
    
    char k_1str[5];
    sprintf(k_1str, "%d", k-1);
    char k_2str[5];
    sprintf(k_2str, "%d", k-2); 

    for(int speciesID=0; speciesID<2; speciesID++)
    {
      //cout << "==== HAO: load the kmer count for (k-1) = " << k_1str << " ==" << endl;
      char kmerFilek_1PathName[1000];
      strcpy(kmerFilek_1PathName, inputFileDir[speciesID]);
      //strcat(kmerFilek_1PathName, "/kmerCount/");
      strcat(kmerFilek_1PathName, speciesName[speciesID]);
      strcat(kmerFilek_1PathName, "_k");
      strcat(kmerFilek_1PathName, k_1str);
      strcat(kmerFilek_1PathName, "_singleStrand_wordcount");
      //cout << "kmerFile: " << kmerFilek_1PathName << endl;
      loadKmerCountHash(kmerFilek_1PathName, speciesID, "kmerCount_1");
      // load the kmer count for MC: k=order+1
      //cout << "==== HAO: load the kmer count for (k-2) = " << k_2str << " ==" << endl;
      char kmerFilek_2PathName[1000];
      strcpy(kmerFilek_2PathName, inputFileDir[speciesID]);
      //strcat(kmerFilek_2PathName, "/kmerCount/");
      strcat(kmerFilek_2PathName, speciesName[speciesID]);
      strcat(kmerFilek_2PathName, "_k");
      strcat(kmerFilek_2PathName, k_2str);
      strcat(kmerFilek_2PathName, "_singleStrand_wordcount");
      //cout << "kmerFile: " << kmerFilek_2PathName << endl;
      loadKmerCountHash(kmerFilek_2PathName, speciesID, "kmerCount_2");
    }

    
    //cout << "==== compute the Hao statistics ==== " << endl;
    double HAOvalue = HAOcompute(ZI, k);
    cout << "hao, " <<  HAOvalue << endl;
    //fout << "hao, " << HAOvalue << endl;
    
  }
  
  
  
  

  // 4. compute JS distance 
  //cout << "== compute the JSdistance == " << endl;
  double entropyRate1;
  double entropyRate2;
  double entropyRateAve;
  double JSdivergence;
  double JSdistance;
  
  if(order[0] != order[1])
  {
    //cout << "ERROR: order 1 MUST equal order 2." << endl;
  }else{
    
    if(order[0] == 0){
      
      vector<SCIENTIFIC_NUMBER> iidProb1;
      vector<SCIENTIFIC_NUMBER> iidProb2;
      vector<SCIENTIFIC_NUMBER> iidProbAve;
      for(int currentIndex = 0; currentIndex < ZI; currentIndex++)
      {
        iidProb1.push_back(HashPw[0][currentIndex]);
        iidProb2.push_back(HashPw[1][currentIndex]);
        iidProbAve.push_back(SciMultiple(SciAddition(HashPw[0][currentIndex], HashPw[1][currentIndex]), 0.5));
      }
      
      entropyRate1 = JScomputeIID(iidProb1);
      entropyRate2 = JScomputeIID(iidProb2);
      entropyRateAve = JScomputeIID(iidProbAve);
      JSdivergence = entropyRateAve - 0.5 * entropyRate1 - 0.5 * entropyRate2;
      JSdistance = pow(JSdivergence, 0.5);
      
      
    }else{
      
      // order1 MUST equal order2
      entropyRate1 = JScomputeMC(iniProb[0], transMatrix[0]);
      entropyRate2 = JScomputeMC(iniProb[1], transMatrix[1]);
      
      vector<SCIENTIFIC_NUMBER> iniProbAve;
      for(int currentIndex = 0; currentIndex < pow(ZI, order[0]); currentIndex++)
      { 
        iniProbAve.push_back(SciMultiple(SciAddition(iniProb[0][currentIndex], iniProb[1][currentIndex]), 0.5));
        //cout << TransToReal(iniProbAve[currentIndex]) << endl;
      }
      
      vector< vector<SCIENTIFIC_NUMBER> > transMatrixAve(pow(ZI,order[0]), vector<SCIENTIFIC_NUMBER>(ZI));
      for(int currentRow = 0; currentRow < pow(ZI, order[0]); currentRow++)
      {
        for(int currentCol = 0; currentCol < ZI; currentCol++)
        {
          // compute JS distance
          transMatrixAve[currentRow][currentCol] = SciMultiple(SciAddition(transMatrix[0][currentRow][currentCol], transMatrix[1][currentRow][currentCol]), 0.5);
          
        }
      }
      entropyRateAve = JScomputeMC(iniProbAve, transMatrixAve);
      
      JSdivergence = entropyRateAve - 0.5 * entropyRate1 - 0.5 * entropyRate2;
      JSdistance = pow(JSdivergence, 0.5);
      
    }
    
    //cout << entropyRate1 << ", " << entropyRate2 << ", ";
    //cout << entropyRateAve << ", " ;
    //cout << JSdivergence << ", " ;
    cout << "JS, " << JSdistance << endl ;
    //fout << "JS, " << JSdistance << endl ;
    
  }
  

  
  // 5. compute Eu, Ma, Ch
  //cout << "== compute the EuMaCh distances == " << endl;
  vector<double> EuMaChDistvalues = EuMaChDistNGS(ZI, k) ;
  cout << "Eu, " << EuMaChDistvalues[0] << endl ;
  cout << "Ma, " << EuMaChDistvalues[1] << endl ;
  cout << "Ch, " << EuMaChDistvalues[2] << endl ;
  //fout << "Eu, " << EuMaChDistvalues[0] << endl ;
  //fout << "Ma, " << EuMaChDistvalues[1] << endl ;
  //fout << "Ch, " << EuMaChDistvalues[2] << endl ;
	
	// 6. compute modified Eu
	//cout << "== compute the modified Eu distances == " << endl;
	vector<double> EuFDistvalues = modifiedEuDistNGS(ZI, k) ;
	cout << "EuF, " << EuFDistvalues[0] << endl;
	//fout << "EuF, " << EuFDistvalues[0] << endl;
  
  // 7. compute Willner
  //cout << "== compute the Willner distance == " << endl;
  if( k > 4)
  {
    cout << "ERROR: no definition of Willner for k>4." << endl;
    
  }else{
    double willner = WillnerDiNGS(ZI, k) ;
    cout << "willner, " << willner << endl ;
    //fout << "willner, " << willner << endl ;
  }
  
  
  
  //fout.close();
  

  return 0;
}

