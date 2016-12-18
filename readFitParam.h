// This small class is used to read fit parameters from text file.
// The keyWord is restricted by user. This function is only for special purpose.

class param{

 public:

  param(string, string);
  float value(string);

 private:

  string flavor_;
  string btag_;

};

param::param(string flavor, string btag){

  flavor_ = flavor;
  btag_   = btag;  

}

float param::value(string keyWord){

  ifstream tFile("/afs/cern.ch/user/v/vieri/work/ZH_Analysis/ZpZHllbb_13TeV/mZhFitParam.txt");

  string tFlavor, tBtag;
  float a_domSb, b_domSb, a_domSg, b_domSg, a_sub1Sb, b_sub1Sb, a_sub1Sg, b_sub1Sg, a_sub2Sb, b_sub2Sb, a_sub2Sg, b_sub2Sg, a_dataSb, b_dataSb;
  float j_data, j_mc, thisParam;

  // ignore first line of text file

  tFile.ignore(1000,'\n');

  while( tFile    >> 
	 tFlavor  >> 
	 tBtag    >> 
	 a_domSb  >> 
	 b_domSb  >> 
	 a_domSg  >> 
	 b_domSg  >> 
	 a_sub1Sb >>
	 b_sub1Sb >> 
	 a_sub1Sg >> 
	 b_sub1Sg >>
	 a_sub2Sb >> 
	 b_sub2Sb >>
	 a_sub2Sg >> 
	 b_sub2Sg >>
	 a_dataSb >> 
	 b_dataSb >> 
	 j_data   >>
	 j_mc     ){
  
    if( flavor_ == tFlavor && btag_ == tBtag ){
    
      if     ( keyWord == "a_domSb"    ) thisParam = a_domSb; 
      else if( keyWord == "a_domSbMin" ) thisParam = a_domSb*0.5;
      else if( keyWord == "a_domSbMax" ) thisParam = a_domSb*1.5;

      else if( keyWord == "b_domSb"    ) thisParam = b_domSb;
      else if( keyWord == "b_domSbMin" ) thisParam = b_domSb*0.5;
      else if( keyWord == "b_domSbMax" ) thisParam = b_domSb*1.5;

      else if( keyWord == "a_domSg"    ) thisParam = a_domSg;
      else if( keyWord == "a_domSgMin" ) thisParam = a_domSg*0.5;
      else if( keyWord == "a_domSgMax" ) thisParam = a_domSg*1.5;

      else if( keyWord == "b_domSg"    ) thisParam = b_domSg;
      else if( keyWord == "b_domSgMin" ) thisParam = b_domSg*0.5;
      else if( keyWord == "b_domSgMax" ) thisParam = b_domSg*1.5;

      else if( keyWord == "a_sub1Sb"    ) thisParam = a_sub1Sb;
      else if( keyWord == "a_sub1SbMin" ) thisParam = a_sub1Sb*0.5;
      else if( keyWord == "a_sub1SbMax" ) thisParam = a_sub1Sb*1.5;

      else if( keyWord == "b_sub1Sb"    ) thisParam = b_sub1Sb;
      else if( keyWord == "b_sub1SbMin" ) thisParam = b_sub1Sb*0.5;
      else if( keyWord == "b_sub1SbMax" ) thisParam = b_sub1Sb*1.5;

      else if( keyWord == "a_sub1Sg"    ) thisParam = a_sub1Sg;
      else if( keyWord == "a_sub1SgMin" ) thisParam = a_sub1Sg*0.5;
      else if( keyWord == "a_sub1SgMax" ) thisParam = a_sub1Sg*1.5;

      else if( keyWord == "b_sub1Sg"    ) thisParam = b_sub1Sg;
      else if( keyWord == "b_sub1SgMin" ) thisParam = b_sub1Sg*0.5;
      else if( keyWord == "b_sub1SgMax" ) thisParam = b_sub1Sg*1.5;

      else if( keyWord == "a_sub2Sb"    ) thisParam = a_sub2Sb;
      else if( keyWord == "a_sub2SbMin" ) thisParam = a_sub2Sb*0.5;
      else if( keyWord == "a_sub2SbMax" ) thisParam = a_sub2Sb*1.5;
      
      else if( keyWord == "b_sub2Sb"    ) thisParam = b_sub2Sb;
      else if( keyWord == "b_sub2SbMin" ) thisParam = b_sub2Sb*0.5;
      else if( keyWord == "b_sub2SbMax" ) thisParam = b_sub2Sb*1.5;

      else if( keyWord == "a_sub2Sg"    ) thisParam = a_sub2Sg;
      else if( keyWord == "a_sub2SgMin" ) thisParam = a_sub2Sg*0.5;
      else if( keyWord == "a_sub2SgMax" ) thisParam = a_sub2Sg*1.5;
      
      else if( keyWord == "b_sub2Sg"    ) thisParam = b_sub2Sg;
      else if( keyWord == "b_sub2SgMin" ) thisParam = b_sub2Sg*0.5;
      else if( keyWord == "b_sub2SgMax" ) thisParam = b_sub2Sg*1.5;

      else if( keyWord == "a_dataSb"    ) thisParam = a_dataSb;
      else if( keyWord == "a_dataSbMin" ) thisParam = a_dataSb*0.5;
      else if( keyWord == "a_dataSbMax" ) thisParam = a_dataSb*1.5;

      else if( keyWord == "b_dataSb"    ) thisParam = b_dataSb;
      else if( keyWord == "b_dataSbMin" ) thisParam = b_dataSb*0.5;
      else if( keyWord == "b_dataSbMax" ) thisParam = b_dataSb*1.5;

      else if( keyWord == "j_data"      ) thisParam = j_data;
      else if( keyWord == "j_dataMin"   ) thisParam = j_data*0.5;
      else if( keyWord == "j_dataMax"   ) thisParam = j_data*1.5;

      else if( keyWord == "j_mc"        ) thisParam = j_mc;
      else if( keyWord == "j_mcMin"     ) thisParam = j_mc*0.5;
      else if( keyWord == "j_mcMax"     ) thisParam = j_mc*1.5;

      else thisParam = 0;
    
    }

  }

  return thisParam;

}
