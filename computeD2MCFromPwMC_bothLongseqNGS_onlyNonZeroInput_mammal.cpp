
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include<iostream>
#include<fstream>
#include<cstring>
//#include<cmath>
#include<vector>
#include <sstream>
#include <string>
#include <stdlib.h> strtol(s.c_str(),0,10);
using namespace std;

# include <stdio.h>
# include <stdlib.h>

#include<algorithm>
#include<ext/hash_map>
using __gnu_cxx::hash_map;
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
hash_map<unsigned long,unsigned long> HashTable[2];
// kmer length = order
hash_map<unsigned long,unsigned long> HashTableOrder[2];
// kmer length = order + 1
hash_map<unsigned long,unsigned long> HashTableOrder_1[2];

//receptacle of Pw(probability of a kmer word) 
hash_map<unsigned long,SCIENTIFIC_NUMBER > HashPw[2];

unsigned long totalKmer[2];
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





void loadPwMCHash(char *currentPwMCFilePathName, int speciesID)
{
  ifstream fin(currentPwMCFilePathName);  
  string currentPwMCLine;
  while(getline(fin, currentPwMCLine, '\n'))
  { 
    unsigned long currentKmerID; double currentPw;
    std::istringstream ss(currentPwMCLine);
    std::string token;
    int colCount = 0;
    while(std::getline(ss, token, ',')) {
      colCount++;
      if(colCount == 1){
        currentKmerID = atoi(token.c_str());
      }else if(colCount == 2){
        currentPw = atof(token.c_str()) ;
        HashPw[speciesID][currentKmerID * ZI + 0] = TransToScientific(currentPw);
      }else if(colCount == 3)
      {
        currentPw = atof(token.c_str()) ;
        HashPw[speciesID][currentKmerID * ZI + 1] = TransToScientific(currentPw);
      }else if(colCount == 4)
      {
        currentPw = atof(token.c_str()) ;
        HashPw[speciesID][currentKmerID * ZI + 2] = TransToScientific(currentPw);
      }else{
        currentPw = atof(token.c_str()) ;
        HashPw[speciesID][currentKmerID * ZI + 3] = TransToScientific(currentPw);
      }
      
      //std::cout << token << '\n';
    }
    
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
  hash_map<unsigned long,SCIENTIFIC_NUMBER > probIID;
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



void pwMC(int ZI, int k, int order, int speciesID)
{
  // compute the transition probability matrix
  int transColSize = ZI;
  int transRowSize = pow(ZI, order);
  vector< vector<SCIENTIFIC_NUMBER> > transMatrix(transRowSize, vector<SCIENTIFIC_NUMBER>(transColSize));
  vector<SCIENTIFIC_NUMBER> iniProb;
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
      transMatrix[currentRow][currentCol] = SciMultiple(probAboveSci,inv_probBelowSci);
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











///////////////////////////////////////////////////////////////////////////////////////////
//KmerCount.out [k] [sample data file] 
int main(int argc, char *argv[])   //EDIT main(int argc, char *argv[])
{
//  HashTable.clear();
//  HashPw.clear();
  //cout << "test" << endl;
  char *speciesName[2];
  vector<int> order(2,0);
  speciesName[0] = argv[1];
  order[0] = atoi(argv[2]);
  speciesName[1] = argv[3];
  order[1] = atoi(argv[4]);
  
  int k = atoi(argv[5]);
  char *kstr = argv[5]; 
  
  int d = atoi(argv[6]);
  char *dstr = argv[6];   
  
  char *dirPath = argv[7];
  
//  char orderstr[5]; 
//  sprintf(orderstr, "%d", order);
//  char order1str[5]; 
//  sprintf(order1str, "%d", (order+1));

  

  // load the kmer count hashtables of the two dataset

  for(int speciesID=0; speciesID<2; speciesID++)
  {
    cout << "======= species " << speciesID << " order " << order[speciesID] << endl;
    char orderstr[5];
    sprintf(orderstr, "%d", order[speciesID]);
    char order1str[5];
    sprintf(order1str, "%d", (order[speciesID]+1));
    
    ///panfs/cmb-panasas2/renj/10species/rn5/rn5_d1_k14_d0_singleStrand_wordcount
    // load the kmer count
    cout << "== load the kmer count ==" << endl;
    char kmerFilePathName[1000];
    strcpy(kmerFilePathName, dirPath);
    strcat(kmerFilePathName, "/");
    strcat(kmerFilePathName, speciesName[speciesID]);
    strcat(kmerFilePathName, "/");
    strcat(kmerFilePathName, speciesName[speciesID]);
    strcat(kmerFilePathName, "_d");
    strcat(kmerFilePathName, dstr);    
    strcat(kmerFilePathName, "_k");
    strcat(kmerFilePathName, kstr);
    strcat(kmerFilePathName, "_d0_singleStrand_wordcount");
    cout << "kmerFile: " << kmerFilePathName << endl;
    loadKmerCountHash(kmerFilePathName, speciesID, "kmerCount");
    
    // load the Pw for kmer under MC model
    cout << "== load the Pw for kmer words, under a MC model with order " << order[speciesID] << " ==" << endl;
    char currentPwMCFilePathName[1000];
    strcpy(currentPwMCFilePathName, dirPath);
    strcat(currentPwMCFilePathName, "/");
    strcat(currentPwMCFilePathName, speciesName[speciesID]);
    strcat(currentPwMCFilePathName, "/");
    strcat(currentPwMCFilePathName, speciesName[speciesID]);
    strcat(currentPwMCFilePathName, "_d");
    strcat(currentPwMCFilePathName, dstr);
    strcat(currentPwMCFilePathName, "_k");
    strcat(currentPwMCFilePathName, kstr);
    strcat(currentPwMCFilePathName, "_order");
    strcat(currentPwMCFilePathName, orderstr);
    strcat(currentPwMCFilePathName, "_singleStrand_pwMC");
    cout << "pwMCFile: " << currentPwMCFilePathName << endl;
    loadPwMCHash(currentPwMCFilePathName, speciesID);
      
  }
    
  
  // compute the D2 statistics
  cout << "== compute the D2 statistics, NOT consider complementary == " << endl;
  vector<double> D2C2Longseqvalues = D2C2computeLongseq(ZI, k);
    
  char outputfile[1000];
  char orderstrA[5];
  sprintf(orderstrA, "%d", order[0]);
  char orderstrB[5];
  sprintf(orderstrB, "%d", order[1]);
  
  strcpy(outputfile, dirPath);
  strcat(outputfile, "/NGS/d");
  strcat(outputfile,dstr);
  strcat(outputfile,"/k");
  strcat(outputfile,kstr);
  strcat(outputfile,"_");
  strcat(outputfile,speciesName[0]);
  strcat(outputfile,"_order");
  strcat(outputfile,orderstrA);
  strcat(outputfile,"_");
  strcat(outputfile,speciesName[1]);
  strcat(outputfile,"_order");
  strcat(outputfile,orderstrB);
  strcat(outputfile,"_d");
  strcat(outputfile,dstr); 
  strcat(outputfile,"_longseq_d2c2score");
  ofstream fout(outputfile);
  for(int i = 0; i < D2C2Longseqvalues.size(); i++)
  {
    //cout << i << endl;
    fout << D2C2Longseqvalues[i] << ", " << endl;
    cout << D2C2Longseqvalues[i] << ", " << endl;
  }
  fout.close();
  
  
  
  cout << "== compute the D2 statistics, consider complementary == " << endl;
  vector<double> D2C2NGSvalues = D2C2computeNGS(ZI, k);
  //printFour(D2C2values);
 
  char outputfile1[1000];
  strcpy(outputfile1, dirPath);
  strcat(outputfile1, "/NGS/d");
  strcat(outputfile1, dstr);
  strcat(outputfile1,"/k");
  strcat(outputfile1,kstr);
  strcat(outputfile1,"_");
  strcat(outputfile1,speciesName[0]);
  strcat(outputfile1,"_order");
  strcat(outputfile1,orderstrA);
  strcat(outputfile1,"_");
  strcat(outputfile1,speciesName[1]);
  strcat(outputfile1,"_order");
  strcat(outputfile1,orderstrB);
  strcat(outputfile1,"_d");
  strcat(outputfile1,dstr); 
  strcat(outputfile1,"_NGS_d2c2score");
  ofstream fout1(outputfile1);
  for(int i = 0; i < D2C2NGSvalues.size(); i++)
  {
    //cout << i << endl;
    fout1 << D2C2NGSvalues[i] << ", " << endl;
    cout << D2C2NGSvalues[i] << ", " << endl;
  }
  fout1.close();
  
  
  

  return 0;
}

