#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include "parser.h"
#ifndef LOOKUP_H
#include "lookup.h"
#define LOOKUP_H
#endif
#include "makeLookup.h"
#include "newmatap.h"
#include "newmatio.h"

#define MIN_READ_DEPTH_INDEL 10
#define N_MAX_ALLELES 20

// parameters for the indel mutation rate model
const float INSERTION_SLOPE = -0.2994;
const float INSERTION_INTERCEPT = -22.8689;
const float DELETION_SLOPE = -0.2856;
const float DELETION_INTERCEPT = -21.9313;



using namespace std;

// Calculate DNM and Null PP
void trio_like_VNTR(indel_t *child,indel_t *mom, indel_t *dad, int flag, 
                     vector<vector<string > > & tgtIndel, 
                     lookup_indel_t & lookupIndel, double mu_scale, 
                     string op_vcf_f, ofstream& fo_vcf, double pp_cutoff, 
                     int RD_cutoff, int& n_site_pass)
{  

  // Read depth filter
  if (child->depth < RD_cutoff ||
      mom->depth < RD_cutoff || dad->depth < RD_cutoff || child->n_alleles > N_MAX_ALLELES) {
    cerr<<"Failed RD filter: "<<child->chr<<":"<<child->pos<<"\t"<<child->depth<<"\t"<<mom->depth<<"\t"<<dad->depth<<endl;
    free(mom->lk);
    free(child->lk);
    free(dad->lk);
    return;
  }

  //cout<<"\nn_alleles: "<<mom->n_alleles<<"\t"<<dad->n_alleles<<"\t"<<child->n_alleles;
  /*if(child->n_alleles > N_MAX_ALLELES) {
    cerr<<"Number of alleles gt max alleles: "<<child->chr<<":"<<child->pos<<endl;   
    return;
  }*/
  //cout<<"0\n";
  n_site_pass += 1;
  Real* a = new Real[sizeof(Real)*mom->gt_len];   
  Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom, sum_denovo, pp_denovo2;       
  Matrix M(1, mom->gt_len);
  Matrix C(child->gt_len,1);
  Matrix D(dad->gt_len,1);
  Matrix P(mom->gt_len, dad->gt_len);
  Matrix F(mom->gt_len*dad->gt_len, child->gt_len);
  Matrix NORMAL(mom->gt_len*dad->gt_len, child->gt_len);
  Matrix DENOVO(mom->gt_len*dad->gt_len, child->gt_len);
  Matrix PP(9,3);
  int i,j,k,l;
  int coor = child->pos;
  char ref_name[50];
  vector<vector<string > > gt_string;
  strcpy(ref_name, child->chr); // Name of the reference sequence 
  //cout<<"1\n";
  cout<<"\nGT len 1: "<<child->gt_len<<endl;
  gt_string.resize(child->gt_len * child->gt_len);
  cout<<"\n2\n";
  //for(int j = 0; j < child->gt_len * child->gt_len; j++) { //dad                                                                                                           //    gt_string[j].resize(child->gt_len);
  //}
  makeVNTRDenovoMask(mom->n_alleles, NORMAL, DENOVO, gt_string); //set the denovo mask, 1 if denovo configuration 0 if null configuration.
  //cout<<"2\n";
  //Load likelihood vectors    
  for (j = 0; j < child->gt_len; ++j) { 
    //for (j = 0; j != 3; ++j) { 
    //cout<<"mom lik "<<mom->lk[j]<<"\t"<<child->gt_len<<endl;
    a[j]=pow(10, -mom->lk[j]/10.); 
  }
  free(mom->lk);
  M<<a;

  for (j = 0; j < child->gt_len; ++j) 
    //for (j = 0; j != 3; ++j) 
    a[j]=pow(10, -dad->lk[j]/10.);
  free(dad->lk);
  D<<a;

  for (j = 0; j < child->gt_len; ++j) 
  //for (j = 0; j != 3; ++j) 
    a[j]=pow(10, -child->lk[j]/10.);
  free(child->lk);
  C<<a;

  delete[] a;

  //cout<<"3\n";

  P = KP(M,D);    
  F = KP(P,C);    
  
  //Find max likelihood of null configuration
  PP = SP(F, NORMAL);   //zeroes out configurations with mendelian error
  maxlike_null = PP.maximum2(i,j);       

  //Find max likelihood of de novo trio configuration
  PP = SP(F, DENOVO);   //zeroes out configurations with mendelian inheritance
  maxlike_denovo = PP.maximum2(k,l); 
  sum_denovo = PP.sum();

  //cout<<"4\n";
  /*cout<<endl<<"i"<<i<<"j "<<j<<"k "<<k<<"l "<<l<<endl;
  cout<<endl<<gt_string.size()<<"\t"<<gt_string[0].size()<<"\t"<<child->gt_len<<endl;

  for(i=0; i<child->gt_len*child->gt_len; i++)
      for(j=0; j<child->gt_len; j++)
      cout<<gt_string[i][j];*/


  /*if(child->gt_len ==6) {
    cout<<"NORMAL"<<endl;
    cout<<NORMAL<<endl;
    }*/
  /*if(child->gt_len == 3) {  
    cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<endl;
    cout<<"NORMAL"<<endl;
    cout<<NORMAL;
    cout<<"DENOVO"<<endl;
    cout<<DENOVO;
    cout<<"M"<<endl;
    cout<<M;
    cout<<"D"<<endl;
    cout<<D;
    cout<<"C"<<endl;
    cout<<C;
    cout<<"P"<<endl;
    cout<<P;
    cout<<"F"<<endl;
    cout<<F;
    cout<<"PP NULL"<<endl;
    cout<<SP(F, NORMAL);
    cout<<"PP DENOVO"<<endl;
    cout<<PP;
    cout<<"ALT "<<child->alt<<endl;
    }*/
  
  //make proper posterior probs
  denom = F.sum();
  pp_denovo = maxlike_denovo/denom; // denovo posterior probability
  pp_denovo2 = sum_denovo/denom; // denovo posterior probability
  pp_null = 1 - pp_denovo; // null posterior probability
  cerr<<"c "<<ref_name<<"\t"<<coor<<endl;
  // Check for PP cutoff
  //if ( pp_denovo > pp_cutoff ) {
  
  //remove ",X" from alt, helps with VCF op.
  string alt = mom->alt;  
  size_t start = alt.find(",X");
  if(start != std::string::npos)
    alt.replace(start, 2, "");
  
  //cout<<"5\n";
  
  cout<<"DENOVO-INDEL child id: "<<child->id;
  cout<<" ref_name: "<<ref_name<<" coor: "<<coor<<" ref_base: "<<mom->ref_base<<" ALT: "<<alt;
  cout<<" F_sum: "<<denom;
  cout<<" maxlike_null: "<<maxlike_null<<" pp_null: "<<pp_null<<" tgt_null: "<<gt_string[i-1][j-1];
  //cout<<" snpcode: "<<lookupIndel.snpcode(i,j)<<" code: "<<lookupIndel.code(i,j);
  cout<<" maxlike_dnm: "<<maxlike_denovo<<" sum_dnm: "<<sum_denovo;
  cout<<" pp_dnm: "<<pp_denovo<<" pp_dnm2: "<<pp_denovo2<<" tgt_denovo: "<<gt_string[k-1][l-1];;
  //cout<<" tgt: "<<tgtIndel[k-1][l-1]<<" lookup: "<<lookupIndel.code(k,l)<<" flag: "<<flag;
  cout<<" READ_DEPTH child: "<<child->depth<<" dad: "<<dad->depth<<" mom: "<<mom->depth;
  cout<<" MAPPING_QUALITY child: "<<child->rms_mapQ<<" dad: "<<dad->rms_mapQ<<" mom: "<<mom->rms_mapQ;
  cout<<endl;
  return;
    /*if(op_vcf_f != "EMPTY") {
      fo_vcf<<ref_name<<"\t";
      fo_vcf<<coor<<"\t";
      fo_vcf<<".\t";// Don't know the rsID
      fo_vcf<<mom->ref_base<<"\t";
      fo_vcf<<alt<<"\t";
      fo_vcf<<"0\t";// Quality of the Call
      fo_vcf<<"PASS\t";// passed the read depth filter
      fo_vcf<<"RD_MOM="<<mom->depth<<";RD_DAD="<<dad->depth;
      fo_vcf<<";MQ_MOM="<<mom->rms_mapQ<<";MQ_DAD="<<dad->rms_mapQ;  
      fo_vcf<<";INDELcode="<<lookupIndel.snpcode(i,j)<<";\t";
      fo_vcf<<"NULL_CONFIG(child/mom/dad):PP_NULL:DNM_CONFIG(child/mom/dad):PP_DNM:RD:MQ\t";
      fo_vcf<<tgtIndel[i-1][j-1]<<":"<<pp_null<<":"; 
      fo_vcf<<tgtIndel[k-1][l-1]<<":"<<pp_denovo<<":"<<child->depth<<":"<<child->rms_mapQ; 
      fo_vcf<<"\n";
      }*/
    //}    

}
