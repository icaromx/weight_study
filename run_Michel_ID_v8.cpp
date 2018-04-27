#include "TROOT.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#define PI 3.14159
#define euler 2.71828

//#include "/grid/fermiapp/products/larsoft/eigen/v3_3_3/include/eigen3/Eigen/Dense"

#include "/usr/local/Cellar/eigen/3.3.4/include/eigen3/Eigen/Dense" //Needed on MACOS
using namespace std;

struct Point {
	float x;
	float y;
	float z;
	float q;
};
struct PCAResults {
	TVector3 centroid;
	pair<TVector3,TVector3> endPoints;
	float length;
	TVector3 eVals;
	vector<TVector3> eVecs;
};
struct TrkPoint{
	double c;
	double x;
	double y;
	double z;
	double q;
};
struct by_y { 
	bool operator()(TrkPoint const &a, TrkPoint const &b) { 
		if(a.y == b.y) return a.x > b.x;
		else return a.y > b.y;
	}
};
struct reverse_by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x < b.x;
        else return a.y < b.y;
    }
};
typedef vector<TrkPoint> track_def;
typedef vector<Point> PointCloud;
void LoadPointCloud(PointCloud &points, const track_def &ord_trk);
PCAResults DoPCA(const PointCloud &points);
double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2);
vector<double> Unit_Vec(double x1,double y1,double z1);
vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2);
double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2);

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////MAIN PROGRAM STARTS////////////////////////////////////////

int main(int argc, char **argv){


	int track_num = atoi(argv[2]);
	///////////////////
	//Define Parameters
	int sample_freq = 40;
	//Ordering Algorithm Parameters
	double alpha = 5.;
	unsigned ang_points = 8;
	int min_points_trk = ang_points*2;
	//double min_cone_ang = 0.75; //41.4 deg
	double min_cone_ang = 0.7; // 49.45 deg
	double back_ang = -0.4;
	double phi = acos(min_cone_ang) * 180. / PI;

	//Cluster Algorithm Parameters
	double EVR = 0.9;
	double min_costheta = -0.85;
	double theta = acos(min_costheta) * 180./ PI;
	double max_costheta =  0.96;
	double prev_win_size = 10;
	double post_win_size = 10;
	double EVR_100 = EVR*100;
	double weight;//, eps_ang, eps_dist;


	std::vector<double> EDIST, EANG;
	int edist_count = 10;
	int eang_count = 10;
	double min_eang = 0.7;//~45 deg
	double max_edist = 4.;
  	for (int i = 0; i <= edist_count; ++i){
  		EDIST.push_back(i * (max_edist / edist_count));
  	}
  	for (int i = 0; i <= eang_count; ++i){
  		EANG.push_back(1 + i * ((min_eang - 1)/eang_count));
  	}

  	///////////////////
  	//Define root output file
  	TFile *f_output;
	if(track_num == -1){
  		f_output = TFile::Open("Weight_Study.root","RECREATE");
	}else{
  		f_output = TFile::Open(Form("results_%s_weight.root",argv[2]),"RECREATE");
	}

	TNtuple *nt_study = new TNtuple("nt_study","nt_study","alpha:Eps_Dist:Eps_Ang:purity:efficiency:ver_rms:ver_mean:t_after_ord_alg:t_sel_as_mich:t_correct_sel_as_mich");
  	/*
  	TNtuple *nt_trk_info = new TNtuple("nt_trk_info","nt_trk_info","run_num:ev_num:cluster_id");
  	TNtuple *nt_results = new TNtuple("nt_results","nt_results","C_ord:alpha:M_ord:Cone_Cosphi:phi:ord_pca_pts:sel:M_sel:evr:min_Costheta:max_theta:purity:efficiency:ver_rms:ver_mean");
  	TNtuple *nt_brk_pts = new TNtuple("nt_brk_pts","nt_brk_pts","run_num:ev_num:cluster_id:alpha:dist:angle:ford_x:ford_y:ford_z:lalpha_x:lalpha_y:lalpha_z");
  	TNtuple *nt_vectors = new TNtuple("nt_vectors","nt_vectors", "run_num:ev_num:cluster_id:angle:pca_x_s:pca_y_s:pca_z_s:pca_x_e:pca_y_e:pca_z_e:lalpha_x_e:lalpha_y_e:lalpha_z_e");
  	TNtuple *nt_pca_points = new TNtuple("nt_pca_points","nt_pca_points","run_num:ev_num:cluster_id:pca_x:pca_y:pca_z");

  	TNtuple *nt_pca_check_evecs = new TNtuple("nt_pca_check_evecs","nt_pca_check_evecs", "run_num:ev_num:cluster_id:ord_x:ord_y:ord_z:eVecs_x:eVecs_y:eVecs_z:pca_vec_x:pca_vec_y:pca_vec_z");
  	TNtuple *nt_pca_check_points = new TNtuple("nt_pca_check_points","nt_pca_check_points","run_num:ev_num:cluster_id:pca_chunk_x:pca_chunk_y:pca_chunk_z:last_pt");
  	TNtuple *nt_trk_pts = new TNtuple("nt_trk_pts", "nt_trk_points","run_num:ev_num:cluster_id:michel:selected:ordered:x:y:z:point_ord:vertex");
  	TNtuple *nt_trk_pts_ord = new TNtuple("nt_trk_pts_ord", "nt_trk_points_ord","run_num:ev_num:cluster_id:ordered:x:y:z:point_ord:vertex");
  	*/
  	///////////////////
  	//OUTPUT COMMENTS
  	cout << "Using only topological cuts with this selection." << endl;

  	///////////////////
  	//READ IN MICHEL LIST
  	ifstream csv_infile("Michel_candidates_vertex_v7.csv");
  	vector<string> TrackData;
  	std::vector<std::vector<double> > Michel_candidates;
  	std::string mline;
  	while (getline(csv_infile, mline,'\n')){
   		TrackData.push_back(mline); //Get each line of the file as a string
  	}
  	int s = TrackData.size();
  	for (unsigned int i=1; i<s; ++i){
   		std::vector<double> v_michel;
    	std::size_t first_comma = TrackData[i].find(",");      // position of the end of the name of each one in the respective string
    	std::size_t second_comma = TrackData[i].find(",", first_comma + 1);
    	std::size_t third_comma = TrackData[i].find(",", second_comma + 1);
    	std::size_t fourth_comma = TrackData[i].find(",", third_comma + 1);
    	std::size_t fifth_comma = TrackData[i].find(",", fourth_comma + 1);
    	double mrun = std::stod(TrackData[i].substr(0,TrackData[i].size()));
    	double meve = std::stod(TrackData[i].substr(first_comma+1,TrackData[i].size()));
    	double mtrk = std::stod(TrackData[i].substr(second_comma+1,TrackData[i].size()));
    	double mvrX = std::stod(TrackData[i].substr(third_comma+1,TrackData[i].size()));
	    double mvrY = std::stod(TrackData[i].substr(fourth_comma+1,TrackData[i].size()));
	    double mvrZ = std::stod(TrackData[i].substr(fifth_comma+1,TrackData[i].size()));
    	//if(mvrX == -1.0 && mvrY == -1.0 && mvrZ == -1.0){
    	//	continue;
	    //}
    	v_michel.push_back(mrun);
    	v_michel.push_back(meve);
    	v_michel.push_back(mtrk);
    	v_michel.push_back(mvrX);
    	v_michel.push_back(mvrY);
    	v_michel.push_back(mvrZ);
    	Michel_candidates.push_back(v_michel);
   	}
   	int mcand_size = Michel_candidates.size();
  	///////////////////////////////////////////////////////////////////////////////////////////////////
  	///////////////////////////////////////////////////////////////////////////////////////////////////

   	for (int eang = 0; eang < EANG.size(); ++eang){
   		for (int edist = 0; edist < EDIST.size(); ++edist){
   			cout << EANG[eang] << ", " << EDIST[edist] << endl;
		   	///////////////////
		   	// COUNTERS
		   	///////////////////
		   	int small_tracks = 0;
		  	int total_num_tracks = 0;
			int tracks_survived_ord_alg = 0;
			int michels_survived_ord_alg = 0;
			int michel_count = 0;
			int track_selected_as_michel = 0;
			int eigenval_cut = 0;
			int ord_alg_cutouts = 0;
			int ord_alg_cutouts_michels = 0;
		  	int sample_counter = 0;
		  	int pca_count = 0;
			
		   	///////////////////
		   	// STUDY OUTPUTS
		   	///////////////////
		   	double purity, efficiency, ver_rms = 0., ver_mean = 0.;
		  	
			
			std::string line;
			std::ifstream ifs(argv[1]);
			

			while(std::getline(ifs, line)){
				//Boolean Parameters
				bool is_pca_check;

			    // ROOT Settings
			    gROOT->Reset();
			    gErrorIgnoreLevel = kError;

			    //Begin reading ROOT file
			    TString filename;
			    filename.Form("%s",line.c_str());    
			    TFile *infile = new TFile(filename);
			    //Extract Event Metadata
			    TTree *Trun = (TTree*)infile->Get("Trun");
			    Int_t run_num;
			    Int_t ev_num;
			    Trun->SetBranchAddress("runNo",&run_num);
			    Trun->SetBranchAddress("eventNo",&ev_num);
			    Trun->GetEntry(0);
			    //cout << "Looking at run " << run_num << " from event " << ev_num << endl;
			    //Extract Coordinate information
			    TTree *T_charge_cluster = (TTree*)infile->Get("T_charge_cluster_nfc"); 
			    Double_t cluster_id;
			    Double_t qx;
			    Double_t qy;
			    Double_t qz;
			    Double_t qc;
			    T_charge_cluster->SetBranchAddress("qx",&qx);
			    T_charge_cluster->SetBranchStatus("qx", kTRUE);
			    T_charge_cluster->SetBranchAddress("qy",&qy);
			    T_charge_cluster->SetBranchStatus("qy", kTRUE);
			    T_charge_cluster->SetBranchAddress("qz",&qz);
			    T_charge_cluster->SetBranchStatus("qz", kTRUE);
			    T_charge_cluster->SetBranchAddress("qc",&qc);
			    T_charge_cluster->SetBranchStatus("qc", kTRUE);
			    T_charge_cluster->SetBranchAddress("cluster_id", &cluster_id);
			    T_charge_cluster->SetBranchStatus("cluster_id", kTRUE);
			    int all_entries = T_charge_cluster->GetEntries();
			    /////////////////////////////////////////////////////////////
			    //Extract Clusters
			    std::vector<Int_t> clusters;
			    Int_t prev_cval;
		    	for (int i = 0; i < all_entries; ++i){
		      		T_charge_cluster -> GetEntry(i);
		    		if (i == 0){
		        		clusters.push_back(cluster_id);
		        		prev_cval = cluster_id;
		      		}else if(prev_cval != cluster_id){
		        		prev_cval = cluster_id;
		        		clusters.push_back(cluster_id);
		      		}
		      	}
		    	//Looking at tracks individually
		    	int num_clusters = clusters.size();
		    	//Loop through individual clusters
		    	std::vector<int> event_kept_trks;
		    	for (int c = 0; c < num_clusters; ++c){
		    		//PARAMETERS
		    		bool is_michel = false, is_selected = false, is_ordered = false;
		    		int cluster = clusters[c];
		    		//cout << "Looking at Run: " << run_num << ", Event: " << ev_num << ", Cluster: " << cluster << endl;
		    		total_num_tracks += 1;
		      		string clus_id = to_string(ev_num) + to_string(cluster);
		      		int int_clus_id  = stoi(clus_id);
		    
		      		//Determine if cluster is a Michel from hand scan
		      		for (int cand = 0; cand < mcand_size; ++cand){
			    			if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
			    				is_michel = true;
			    				break;
			    			}
			    	}
		    
		      		//Load track info
		      		track_def trk;
		      		//Load every point (cluster id, x, y, z, charge) of a track into the trk object.
		      		for (int i = 0; i < all_entries; ++i){
			        	T_charge_cluster -> GetEntry(i);
			        	//Will only store information for current cluster
			        	if(cluster_id != cluster) continue;
			        	if(track_num != -1){
			        		if(int_clus_id != track_num) continue;
			        	}
				        TrkPoint tempPoint;
				        tempPoint.c = cluster_id;
				        tempPoint.x = qx;
				        tempPoint.y = qy;
				        tempPoint.z = qz;
				        tempPoint.q = qc;
				        trk.push_back(tempPoint);
			      	}
			      	//track size has to be larger than the moving window size
			      	if(trk.size() < min_points_trk + 1){
			      		small_tracks += 1;
			      		continue; // #CUT
			      	}
			      	//////////////////////////////////////////////////
			      	////////ORDERING ALGORITHM BEGINS/////////////////
			      	//Sort track in descending y value
			      	std::sort(trk.begin(), trk.end(), by_y());
			      	int trk_size = trk.size();
			      	track_def points_left;
			      	track_def points_gd;
			      	track_def ord_trk;
			      	//Store points being tested by ordering algorithm
			      	for (int i = 1; i < trk_size; ++i){
			        	TrkPoint tempPoint;
			        	tempPoint.c = trk[i].c;
			        	tempPoint.x = trk[i].x;
			        	tempPoint.y = trk[i].y;
			        	tempPoint.z = trk[i].z;
			        	tempPoint.q = trk[i].q;
			        	points_left.push_back(tempPoint);
			      	}
			      	int pl_size = points_left.size();
			      	std::sort(points_left.begin(), points_left.end(), by_y());
			      	//Highest y value point is the first point in the oredered track
			      	ord_trk.push_back(trk[0]);
				   
				    //Ordering algorithm parameters
				    double old_dist = 10000000.;
				    int low_dist_at = -1;
				    double dist;
				    double low_old_dist = 1000.;
				    double low_ord_y = 2000.;

				    bool cone_test_fail = true;
					bool ongoing_cone_test = true;
					bool closest_point_found_cone = false;
					bool closest_point_found_volume = false;
					bool closest_point_start_found = false;
					bool print_vals = false;

				    int m = 0;
				    int i = 0;
				    int pca_counter = 0;
				    
				    double ang_eigV_x=0., ang_eigV_y=0., ang_eigV_z=0.;
				    double ob_temp_x, ob_temp_y, ob_temp_z;
				    double ftest_point_x, ftest_point_y, ftest_point_z; 
				    double pca_eVec_x, pca_eVec_y, pca_eVec_z; 
		          	double ford_point_x, ford_point_y, ford_point_z, last_dist;
				    double ang_dotP, fangle;
				    std::vector<double> vec_pca_ang;
				    std::vector<double> vec_pca;
				    std::vector<double> vec_ob_temp;
				    std::vector<double> vec_temp;
				    std::vector<double> vec_last_ord;
				    std::vector<double> vec_pca_test;
				    //double pca_end_x, pca_end_y, pca_end_z;
				    double pca_start_x, pca_start_y, pca_start_z;
				    double pca_dotP;
				    std::vector<double> pca_vec;
				    std::vector<double> cone_vec;
				    double cone_dotP, dotP_low_dist;
				    double flip;
				    double ang_decide;
				    track_def pca_points;
				    TrkPoint VertexPoint;


					float xVal;
					float yVal;
		  			float zVal;
		  			track_def ang_chunk;
		  			
			    	while(pl_size != 0){
			    		ang_chunk.clear();
			    		vec_pca_ang.clear();
			    		low_dist_at = -1;
			    		//vec_pca_test.clear();
			    		xVal = 0.;
			    		yVal = 0.;
			    		zVal = 0.;

			    		if(ord_trk.size() > ang_points){
			    			for (unsigned p = ord_trk.size() - ang_points; p < ord_trk.size(); ++p){
				     			TrkPoint ang_tempPoint;
						        ang_tempPoint.c = ord_trk[p].c;
						        ang_tempPoint.x = ord_trk[p].x;
						        ang_tempPoint.y = ord_trk[p].y;
						        ang_tempPoint.z = ord_trk[p].z;
						        ang_tempPoint.q = ord_trk[p].q;
						        ang_chunk.push_back(ang_tempPoint);	
			    			}
				    		PointCloud ang_pointcloud;
				    		PCAResults ang_results;
				    		LoadPointCloud(ang_pointcloud, ang_chunk);
						    ang_results = DoPCA(ang_pointcloud);
							
						    vec_pca_ang.push_back(ang_results.eVecs[0](0));
							vec_pca_ang.push_back(ang_results.eVecs[0](1));
							vec_pca_ang.push_back(ang_results.eVecs[0](2));

							xVal = ang_results.endPoints.first(0);
							yVal = ang_results.endPoints.first(1);
		  					zVal = ang_results.endPoints.first(2);

		  					/*
		  					pca_counter = 0;
			  				while ((((ang_results.endPoints.first(0) < ang_results.endPoints.second(0)) && (xVal < ang_results.endPoints.second(0))) || ((ang_results.endPoints.first(0) >= ang_results.endPoints.second(0)) && (xVal > ang_results.endPoints.second(0)))) && (((ang_results.endPoints.first(1) < ang_results.endPoints.second(1)) && (yVal < ang_results.endPoints.second(1))) || ((ang_results.endPoints.first(1) >= ang_results.endPoints.second(1)) && (yVal > ang_results.endPoints.second(1)))) && (((ang_results.endPoints.first(2) < ang_results.endPoints.second(2)) && (zVal < ang_results.endPoints.second(2))) || ((ang_results.endPoints.first(2) >= ang_results.endPoints.second(2)) && (zVal > ang_results.endPoints.second(2))))) {
		    					if(pca_counter == 0){
		    						pca_start_x = xVal;
		    						pca_start_y = yVal;
		    						pca_start_z = zVal;
		    					}
		    					pca_counter += 1;
		    					xVal += 0.5*ang_results.eVecs[0](0);
		    					yVal += 0.5*ang_results.eVecs[0](1);
		    					zVal += 0.5*ang_results.eVecs[0](2);
		  					}
		  					pca_end_x = xVal;
		  					pca_end_y = yVal;
		  					pca_end_z = zVal;
		  					pca_vec = Unit_Vec_NO(pca_start_x,pca_end_x, pca_start_y, pca_end_y, pca_start_z, pca_end_z);
		  					*/
		  					
		  					pca_vec = Unit_Vec_NO(ang_chunk.at(0).x,ang_chunk.back().x,ang_chunk.at(0).y,ang_chunk.back().y,ang_chunk.at(0).z,ang_chunk.back().z);
		  					pca_dotP = dotProdFunc(pca_vec[0],vec_pca_ang[0],pca_vec[1],vec_pca_ang[1],pca_vec[2],vec_pca_ang[2]);
		  					if (pca_dotP > -pca_dotP){
		  						flip = 1.;
		  					}else{
		  						flip = -1.;
		  					}
		  					vec_pca_ang[0] = flip * vec_pca_ang[0];
		  					vec_pca_ang[1] = flip * vec_pca_ang[1];
		  					vec_pca_ang[2] = flip * vec_pca_ang[2];
					    }else{
					    	pca_vec.push_back(0);
					    	pca_vec.push_back(0);
					    	pca_vec.push_back(0);
					    	vec_pca_ang.push_back(0);
					    	vec_pca_ang.push_back(0);
					    	vec_pca_ang.push_back(0);
					    }

					    cone_test_fail = true;
					    ongoing_cone_test = true;
						closest_point_found_cone = false;
						closest_point_found_volume = false;
						closest_point_start_found = false;
						
						
						dotP_low_dist = 0;

						double old_weight = 0;
						if(ord_trk.size() > ang_points){    
				        	for (int j = 0; j < pl_size; ++j){
				        		cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
				        		cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
				        		if(cone_dotP > min_cone_ang){
				        			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
				        			if(dist < alpha){
						        		weight = pow(euler, -dist/EDIST[edist]) * pow(euler, -min_cone_ang/EANG[eang]);
						        		if(weight > old_weight){
						        			old_weight = weight;
						        			dotP_low_dist = cone_dotP;
					            			low_dist_at = j;
					            			closest_point_found_cone = true;
										}
									}			        	
				        		}
				        	}
				        }
			        	if (closest_point_found_cone == false && ord_trk.size() > ang_points){
			        		for (int j = 0; j < pl_size; ++j){
			        			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
				        		cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
				        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
				        		weight = pow(euler, -dist/EDIST[edist]) * pow(euler, -min_cone_ang/EANG[eang]);
					          	if(dist < alpha){
					          		if(cone_dotP > back_ang){	
						          		if(weight > old_weight){
						            		old_weight = weight;
						            		low_dist_at = j;
						            		closest_point_found_volume = true;
						          		}
						          	}
					          	}	
		        			}
			        	}
			        	if(closest_point_found_cone == false && closest_point_found_volume == false && ord_trk.size() > ang_points){
			        		for (int j = 0; j < pl_size; ++j){
				        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
				          		if (dist < old_dist){
				            		old_dist = dist;
				            		low_dist_at = j;
				          		}
		        			}
			        	}
			        	if(ord_trk.size() <= ang_points){
			        		for (int j = 0; j < pl_size; ++j){
			        			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
				        		cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
				        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
				        		weight = pow(euler, -dist/EDIST[edist]) * pow(euler, -min_cone_ang/EANG[eang]);
				          		if(weight > old_weight){
				            		old_weight = weight;
				            		low_dist_at = j;
				            		closest_point_start_found = true;
				          		}
		        			}
			        	}
			        	//Point with the shortest distance to the last ordered point
			        	TrkPoint tempPoint;
			        	tempPoint.c = points_left[low_dist_at].c;
			        	tempPoint.x = points_left[low_dist_at].x;
			        	tempPoint.y = points_left[low_dist_at].y;
			        	tempPoint.z = points_left[low_dist_at].z;
			        	tempPoint.q = points_left[low_dist_at].q;

			        	ob_temp_x = tempPoint.x - ord_trk.back().x;
			        	ob_temp_y = tempPoint.y - ord_trk.back().y;
			        	ob_temp_z = tempPoint.z - ord_trk.back().z;
			        	vec_ob_temp = Unit_Vec(ob_temp_x,ob_temp_y,ob_temp_z);

			        	if (closest_point_start_found){
			        		ord_trk.push_back(tempPoint);
			          		if (tempPoint.y < low_ord_y){
			          			low_ord_y = tempPoint.y;
			          		}
			          		old_dist = 10000000;
			          		points_left.erase(points_left.begin() + low_dist_at);
			          		pl_size = points_left.size();
			          	}else if(closest_point_found_cone){
			        		//Keep the next point
			        		ord_trk.push_back(tempPoint);
			          		if (tempPoint.y < low_ord_y){
			          			low_ord_y = tempPoint.y;
			          		}
			          		old_dist = 10000000;
			          		points_left.erase(points_left.begin() + low_dist_at);
			          		pl_size = points_left.size();
			          		i = 0;
			          	}else if(closest_point_found_volume){
			        		//Keep the next point
			        		ord_trk.push_back(tempPoint);
			          		if (tempPoint.y < low_ord_y){
			          			low_ord_y = tempPoint.y;
			          		}

			          		old_dist = 10000000;
			          		points_left.erase(points_left.begin() + low_dist_at);
			          		pl_size = points_left.size();
			          		i = 0;
			        	}else{
			        		if(old_dist < low_old_dist){
					    		vec_last_ord.clear();
					    		vec_last_ord.push_back(ord_trk.back().x);
					    		vec_last_ord.push_back(ord_trk.back().y);
					    		vec_last_ord.push_back(ord_trk.back().z);
				                last_dist = old_dist;
				                ftest_point_x = tempPoint.x;
				                ftest_point_y = tempPoint.y;
				                ftest_point_z = tempPoint.z;

				                vec_temp = vec_ob_temp;
				                pca_points = ang_chunk;
				                vec_pca.clear();
				                vec_pca.push_back(vec_pca_ang[0]);
				                vec_pca.push_back(vec_pca_ang[1]);
				                vec_pca.push_back(vec_pca_ang[2]);

				                fangle = cone_dotP;
				                low_old_dist = old_dist;
		              		}
					    	points_gd.push_back(tempPoint);
					        old_dist = 10000000;
					        points_left.erase (points_left.begin() + low_dist_at);
					        pl_size = points_left.size();
					        i++;
			        	}
			        	if (pl_size == 0) break;
			    	}
					double bottom_dist;
					bottom_dist = abs(trk.back().y - low_ord_y);
					// If distance between lowest y value of unordered track and 
					// lowest y value of ordered track is greater than 10 cm.
					fangle = acos(fangle) * 180. / PI;
					//cout << vec_last_ord.size() << endl;
					
					if(bottom_dist > 10.){
			    		if(is_michel) ord_alg_cutouts_michels += 1;
					}else{
						is_ordered = true;
						if(is_michel) michels_survived_ord_alg += 1;
					}
					if(is_ordered) tracks_survived_ord_alg += 1;

					//Finished Ordering Points
				    //////////////////////////

			    	//
					if(is_ordered){
						//////////////////////////
						//Starting Moving Window
					 	track_def prev_chunk;
					 	track_def post_chunk;
					 	PointCloud prev_points;
					 	PointCloud post_points;
					 	//TrkPoint VertexPoint;
					 	//cout << "#################################################################" << endl;
						//cout << "Looking at track " << cluster << " from run = " << run_num << ", event = " << ev_num << endl;
						//cout << "Before window size: " << prev_win_size << "; After window size: " << post_win_size << endl;
					 	double dotProd, min_ang = -360.;
					 	double ev_lowest = 100000;
						double min_prod = 1000000000;
						double prev_eigenratio, post_eigenratio;
						double vertex_prev_ratio, vertex_post_ratio;
					 	track_def kept_prev_chunk;
					 	track_def kept_post_chunk;
					 	double prev_eVecs[3];
					 	double post_eVecs[3];
					 	int vertex;
					 	double prev_angle;
					 	double vertex_res;
					 	unsigned trk_points = 100;
					 	double prev_dotP_decide, post_dotP_decide;
					 	std::vector<double> prev_win_vec, post_win_vec;
					    if(ord_trk.size() < min_points_trk*2) continue;
					    for (int i = prev_win_size + 1; i < ord_trk.size() - post_win_size - 1; ++i){
							track_def prev_chunk;
							track_def post_chunk;
							PointCloud prev_points;
							PointCloud post_points;    
						 	PCAResults prev_results;
						 	PCAResults post_results;
						 	TrkPoint prev_first_point;
						 	TrkPoint post_last_point;
						 	prev_first_point = ord_trk[i - prev_win_size];
						 	post_last_point = ord_trk[i + post_win_size];
							
				     		for (int j = i - prev_win_size; j < i; ++j){
				     			TrkPoint prev_tempPoint;
						        prev_tempPoint.c = ord_trk[j].c;
						        prev_tempPoint.x = ord_trk[j].x;
						        prev_tempPoint.y = ord_trk[j].y;
						        prev_tempPoint.z = ord_trk[j].z;
						        prev_tempPoint.q = ord_trk[j].q;
						        prev_chunk.push_back(prev_tempPoint);	
				     		}

					     	for (int j = i + 1; j < i + post_win_size + 1; ++j){
				     			TrkPoint post_tempPoint;
						        post_tempPoint.c = ord_trk[j].c;
						        post_tempPoint.x = ord_trk[j].x;
						        post_tempPoint.y = ord_trk[j].y;
						        post_tempPoint.z = ord_trk[j].z;
						        post_tempPoint.q = ord_trk[j].q;
						        post_chunk.push_back(post_tempPoint);
					     	}

					     	LoadPointCloud(prev_points, prev_chunk);
					     	LoadPointCloud(post_points, post_chunk);
				     		prev_results = DoPCA(prev_points);
				     		post_results = DoPCA(post_points);
				     		double new_prev_decide, new_post_decide;

				     		prev_win_vec = Unit_Vec_NO(prev_chunk.back().x,prev_chunk.at(0).x,prev_chunk.back().y,prev_chunk.at(0).y,prev_chunk.back().z,prev_chunk.at(0).z);
				     		prev_dotP_decide = dotProdFunc(prev_results.eVecs[0](0),prev_win_vec[0],prev_results.eVecs[0](1),prev_win_vec[1],prev_results.eVecs[0](2),prev_win_vec[2]);		     		
				     		if(prev_dotP_decide > -prev_dotP_decide){
				     			new_prev_decide = -1.;
				     		}else{
				     			new_prev_decide = 1.;
				     		}
				     		post_win_vec = Unit_Vec_NO(post_chunk.back().x,post_chunk.at(0).x,post_chunk.back().y,post_chunk.at(0).y,post_chunk.back().z,post_chunk.at(0).z);
				     		post_dotP_decide = dotProdFunc(post_results.eVecs[0](0),post_win_vec[0],post_results.eVecs[0](1),post_win_vec[1],post_results.eVecs[0](2),post_win_vec[2]);
				     		if(post_dotP_decide > -post_dotP_decide){
				     			new_post_decide = 1.;
				     		}else{
				     			new_post_decide = -1.;
				     		}

				     		prev_results.eVecs[0] = new_prev_decide*prev_results.eVecs[0];
							post_results.eVecs[0] = new_post_decide*post_results.eVecs[0];

							prev_eigenratio = prev_results.eVals(0)/(prev_results.eVals(0) + prev_results.eVals(1) + prev_results.eVals(2));
							post_eigenratio = post_results.eVals(0)/(post_results.eVals(0) + post_results.eVals(1) + post_results.eVals(2));
							dotProd = dotProdFunc(prev_results.eVecs[0](0),post_results.eVecs[0](0),prev_results.eVecs[0](1),post_results.eVecs[0](1),prev_results.eVecs[0](2),post_results.eVecs[0](2));
				     		//dotProd = prev_results.eVecs[0](0)*post_results.eVecs[0](0) + prev_results.eVecs[0](1)*post_results.eVecs[0](1) + prev_results.eVecs[0](2)*post_results.eVecs[0](2);
				     		if(i == prev_win_size){ 
				     			prev_angle = dotProd;
				     		}else{
				     			if(abs(prev_angle - dotProd) > 1.9){
				     				cout << "here" << endl;
				     				prev_angle = dotProd;
				     				continue;
				     			}
				     		}
							if(prev_eigenratio < EVR) continue;
							if(post_eigenratio < EVR) continue;
							if(dotProd > min_ang){
								//cout <<  acos(dotProd) * 180./PI << endl;
								min_ang = dotProd;	
								vertex_prev_ratio = prev_eigenratio;
								vertex_post_ratio = post_eigenratio;
				     			vertex = i;
								VertexPoint.c = ord_trk[vertex].c;
				     			VertexPoint.x = ord_trk[vertex].x;
				     			VertexPoint.y = ord_trk[vertex].y;
				     			VertexPoint.z = ord_trk[vertex].z;
				     			VertexPoint.q = ord_trk[vertex].c;
				     			kept_prev_chunk = prev_chunk;
				     			kept_post_chunk = post_chunk;
				     			prev_eVecs[0] = prev_results.eVecs[0](0);
								prev_eVecs[1] = prev_results.eVecs[0](1);
								prev_eVecs[2] = prev_results.eVecs[0](2);
								post_eVecs[0] = post_results.eVecs[0](0);
								post_eVecs[1] = post_results.eVecs[0](1);
								post_eVecs[2] = post_results.eVecs[0](2);
			     			}
				    	}
				   		//////////////////////////////////////////////////////////////////////
					    double dist_low_vert_y;
					    dist_low_vert_y = abs(VertexPoint.y - low_ord_y);
					    //cout << min_ang << ", " << acos(min_ang) * 180./PI << endl;
					    //Michel electron cutoff distance based on energy spectrum
					    if(is_ordered != true) continue;
					    if(dist_low_vert_y > 15) continue;
					    if(min_ang < min_costheta) continue;
					    if(min_ang > max_costheta) continue;
					    track_selected_as_michel += 1;
					    is_selected = true;
						if(is_michel) michel_count += 1;
				    	for (int cand = 0; cand < mcand_size; ++cand){
				    		if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
				    			if(Michel_candidates[cand][3] == 0 && Michel_candidates[cand][4] == 0 && Michel_candidates[cand][5] == 0) break;
				    			vertex_res = Pythagoras(Michel_candidates[cand][3],VertexPoint.x,Michel_candidates[cand][4],VertexPoint.y,Michel_candidates[cand][5],VertexPoint.z);
				    			ver_rms += pow(vertex_res, 2.);
				    			ver_mean += vertex_res;		    			
				    			break;
				    		}
				    	}
				    }
		  		}
		   		infile->Close();
			}
		    purity = ((float)michel_count)/((float)track_selected_as_michel);
		    efficiency = ((float)michel_count)/((float)(s-1));
		    ver_rms = sqrt((1./((float)michel_count)) * ver_rms);
		    ver_mean = (1./((float)michel_count)) * ver_mean;

		    cout << "###################################################################" << endl;
		    cout << "Total number of tracks = " << total_num_tracks << "; Number of Michels in sample = " << s - 1 << endl;
		    cout << "Tracks smaller than the pca window = " << small_tracks << endl;
		    cout << "###################################################################" << endl;
		    cout << "################### Ordering of Tracks ############################" << endl;
		    cout << "###################################################################" << endl;
		    cout << "Tracks after ordering algorithm = " << tracks_survived_ord_alg << " with alpha = " << alpha << endl;
		    cout << "Michel clusters that survived the ordering algorithm = " << michels_survived_ord_alg << endl;
		    cout << "Michel clusters that were cut out by ordering algorithm = " << ord_alg_cutouts_michels << endl;
		    cout << "Cone Angle Criterion: Cos(/phi) > " << min_cone_ang << " => /phi < " << phi << endl;
		    cout << "Points considered for PCA calculations: " << ang_points << endl;
		    cout << "###################################################################" << endl;
		    cout << "###################### Michel ID part #############################" << endl;
		    cout << "###################################################################" << endl;
		    cout << "Tracks that were selected as Michels  = " << track_selected_as_michel << endl;
		    cout << "Tracks correctly selected as selected as Michels = " << michel_count << endl;
		    cout << "Cuts applied on bk, ak: " << EVR << endl;
		    cout << "Cuts applied on Cos(/theta): Min = " << min_costheta << ", Max = " << max_costheta << endl;
		    cout << "Cuts applied on /theta: Max = " << acos(min_costheta) * 180./ PI << ", Min = " << acos(max_costheta) * 180./ PI << endl;
		    cout << "###################################################################" << endl;
		    cout << "########################## Results #################################" << endl;
		    cout << "###################################################################" << endl;
		    cout << "Purity = " << purity << ", Efficiency = " << efficiency << ", Vertex RMS = " << ver_rms << ", Vertex Mean = " << ver_mean << endl;
		    cout << "####################################################################" << endl;

		    nt_study -> Fill(alpha,EDIST[edist],EANG[eang],purity,efficiency,ver_rms,ver_mean,tracks_survived_ord_alg,track_selected_as_michel,michel_count);
		}
	}
    

    //nt_results->Fill(tracks_survived_ord_alg,alpha,michels_survived_ord_alg,min_cone_ang,phi,ang_points,track_selected_as_michel,michel_count,EVR,min_costheta,theta,purity,efficiency,ver_rms,ver_mean);
      	
  	f_output->Write();
  	f_output->Close();
  	cout << "DONE" << endl;
  	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////END OF MAIN PROGRAM/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////START OF FUNCTIONS//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2){
	double dist;
	dist = sqrt(pow(x2-x1,2.) + pow(y2-y1,2.) + pow(z2-z1,2.));
	return dist;
}

vector<double> Unit_Vec(double x1,double y1,double z1){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,0.0,y1,0.0,z1,0.0);
	v.push_back(x1/norm);
	v.push_back(y1/norm);
	v.push_back(z1/norm);
	return v;
}

vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,x2,y1,y2,z1,z2);
	v.push_back((x2-x1)/norm);
	v.push_back((y2-y1)/norm);
	v.push_back((z2-z1)/norm);
	return v;
}

double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2){
	double dotP;
	dotP = x1*x2 + y1*y2 + z1*z2;
	return dotP;
}


void LoadPointCloud(PointCloud &points, const track_def &ord_trk) {
  for (int i = 0; i < ord_trk.size(); ++i){
    Point tempPoint;
    tempPoint.x = ord_trk.at(i).x;
    tempPoint.y = ord_trk.at(i).y;
    tempPoint.z = ord_trk.at(i).z;
    tempPoint.q = ord_trk.at(i).q;
    points.push_back(tempPoint);

  }
  return;
}

PCAResults DoPCA(const PointCloud &points) {
  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;
  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;
  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }
  if (nThreeDHits == 0) {
    PCAResults results;
    return results; 
  }
  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);
      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);


    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }
  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;
  return results;
}

