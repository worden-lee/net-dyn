/*  
 *  Example usage of igraph package to analyse (Pajek) networks in C code.
 *  
 *  It shows:
 *  - How to read in a network in a foreign format (here: Pajek)
 *  - How to use a vertex selector and a vertex iterator
 *  - How to obtain structural network properties (here: closeness and betweenness centrality)
 *  - How to de-allocate objects
 *
 *  Input: Pajek network file (.net)
 *  Output: Network properties printed to stdout
 *
 *  Example by 					Walter de Back, Collegium Budapest, Hungary

 *  Modified by Lee Worden to input an edgelist rather than Pajek file, and to output to a csv file.


 *	COMPILE: 					gcc analyse-net.c -I/usr/local/igraph -L/usr/local/lib -ligraph -o analyse-net (may need to change paths)


	igraph library. 
	Copyright (C) 2003-2008 Gábor Csárdi and Tamás Nepusz <csardi@rmki.kfki.hu> 
	MTA RMKI, Konkoly-Thege Miklós st. 29-33., Budapest 1121, Hungary.

	DOWNLOAD IGRAPH:			http://cneurocvs.rmki.kfki.hu/igraph/download.html
    IGRAPH REFERENCE MANUAL: 	http://cneurocvs.rmki.kfki.hu/igraph/doc/html/index.html
 
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "/usr/local/include/igraph/igraph.h" // change to path of igraph.h
#include "igraph/igraph.h" // change to path of igraph.h

int main(int argc,char* argv[]){
	
	char infile[FILENAME_MAX], outfile[FILENAME_MAX]; 

	if(argc==3){
		strcpy(infile,  argv[1]);
		strcpy(outfile, argv[2]);
	}else{
		printf("\nUSAGE: ./analyse-net [infile.edgelist] [outfile.csv]\n\n");
		exit(0);
	}

	// Set attribute handler (to allow us to get vertex IDs later on)
	igraph_i_set_attribute_table(&igraph_cattribute_table); // Attribute handler should be set before reading network!

	int i;
	igraph_t graph;

	// Read in Pajek .net file
	FILE* graphfile = fopen(infile,"r");
	int info = -1;
	if(graphfile!=NULL){
	  info = igraph_read_graph_edgelist(&graph, graphfile, 0, IGRAPH_DIRECTED);
	  if(info==0)		
	    printf("file loaded.\n");
	  else{
	    printf("ERROR occured during loading file!\n");
	    exit(-1);
	  }
	}
	fclose(graphfile);

	// Make a vertex selector (vs) that holds all vertices in the graph 
	igraph_vs_t all_vertices;	
	igraph_vs_all(&all_vertices); // all_vertices will contain all vertices of the graph -- NOTE: you can also use igraph_vss_all() 

	// Check how many vertices are in the vertex selector
	igraph_integer_t num_vertices;
	igraph_vs_size(&graph, &all_vertices, &num_vertices);
	printf("\nNumber of Vertices = %d \n",(int)num_vertices);

	// Get names of vertices ('id' in pajek jargon) in an array (shown here with and without vertex iterator) 
	//char* ids[(int)num_vertices]; 

	/*
		Get all IDs of vertices in array of strings using the Vertex iterator

	printf("\nVERTEX NAMES: \n");
	igraph_vit_t vit; // vertex iterator
	igraph_vit_create(&graph, all_vertices, &vit); // initialize the vertex iterator
	while (!IGRAPH_VIT_END(vit)) {
		int p = (int) IGRAPH_VIT_GET(vit);
		ids[p] = igraph_cattribute_VAS(&graph, "id", p);
		printf("Vertex %d =\t%s\n",p,ids[p]);
		IGRAPH_VIT_NEXT(vit);
	}
	*/

	/*
		Do the same (put all vertex-IDs in an array ids[]),
		 but now without the Vertex iterator (which is a bit more readable IMHO)
	*/
/* 	printf("\nVERTEX NAMES: \n");		 */
/* 	for(i=0;i<(int)num_vertices;i++){ */
/* 		ids[i] = (char*) igraph_cattribute_VAS(&graph, "id", i); */
/* 		printf("Vertex %d =\t%s\n",i,ids[i]); */
/* 	} */

	// CENTRALITY (IN)	

	/* Initialize a vector that will hold the results of the network analysis */
	igraph_vector_t in_closeness;
 	igraph_vector_init(&in_closeness, 0); // initialize result vector with length 0 

	info = igraph_closeness(&graph, &in_closeness, igraph_vss_all(), IGRAPH_IN);
	if(info!=0)		
	{ puts("ERROR occured during in closeness!");
	  exit(-1);
	}

	/*
	printf("\nCLOSENESS CENTRALITY (IN) \n");
	for(i=0;i<igraph_vector_size(&in_closeness);i++){
		printf("%d %10f\n",i,igraph_vector_e(&in_closeness,i));
	}
	*/

	// CENTRALITY (OUT)	
	igraph_vector_t out_closeness;
 	igraph_vector_init(&out_closeness, 0); // initialize result vector with length 0 
	info = igraph_closeness(&graph, &out_closeness, igraph_vss_all(), IGRAPH_OUT);
	if(info!=0)		
	{ puts("ERROR occured during out closeness!");
	  exit(-1);
	}

	/*
	printf("\nCLOSENESS CENTRALITY (OUT) \n");		
	for(i=0;i<igraph_vector_size(&out_closeness);i++){
		printf("%d %10f\n",i,igraph_vector_e(&out_closeness,i));		
	}
	*/

	// BETWEENNESS
	igraph_vector_t betweenness;
 	igraph_vector_init(&betweenness, 0); // initialize result vector with length 0 
	info = igraph_betweenness(&graph, &betweenness, igraph_vss_all(), IGRAPH_DIRECTED);
	if(info!=0)		
	{ puts("ERROR occured during betweenness!");
	  exit(-1);
	}

	/*
	printf("\nBETWEENNESS CENTRALITY \n");		
	for(i=0;i<igraph_vector_size(&betweenness);i++){
		printf("%d %10f\n",i,igraph_vector_e(&betweenness,i));		
	}
	*/

	// UNDIRECTED EIGENVECTOR CENTRALITY
	igraph_vector_t und_eigen;
 	igraph_vector_init(&und_eigen, 0); // initialize result vector with length 0 
	igraph_arpack_options_t options;
	igraph_arpack_options_init(&options);
	info = igraph_eigenvector_centrality(&graph, &und_eigen, NULL, 0, NULL, &options);
	if(info!=0)		
	{ puts("ERROR occured during eigenvector centrality!");
	  exit(-1);
	}

#define absolute(x) ((x) < 0 ? -(x):(x))
	/*
	printf("\nUNDIRECTED EIGENVECTOR CENTRALITY \n");		
	for(i=0;i<igraph_vector_size(&und_eigen);i++){
		printf("%d %10f\n",i,absolute(igraph_vector_e(&und_eigen,i)));		
	}
	*/

	// PAGERANK (similar to directed eigenvalue centrality)
	igraph_vector_t pagerank;
	igraph_vector_init(&pagerank,0);
	igraph_real_t eigenvalue;
	igraph_arpack_options_init(&options);
	info = igraph_pagerank(&graph, &pagerank, &eigenvalue, igraph_vss_all(), IGRAPH_DIRECTED,
			       0.85, NULL, &options);
	if(info!=0)		
	{ puts("ERROR occured during pagerank!");
	  exit(-1);
	}
	/*
	printf("\nPAGERANK \n");		
	for(i=0;i<igraph_vector_size(&pagerank);i++){
		printf("%d %10f\n",i,igraph_vector_e(&pagerank,i));		
	}
	*/

	// KLEINBERG HUB SCORE
	igraph_vector_t hub_score;
	igraph_vector_init(&hub_score,0);
	igraph_arpack_options_init(&options);
	info = igraph_hub_score(&graph, &hub_score, NULL, 0, &options);
	if(info!=0)		
	{ puts("ERROR occured during hub score!");
	  exit(-1);
	}

	// KLEINBERG AUTHORITY SCORE
	igraph_vector_t authority_score;
	igraph_vector_init(&authority_score,0);
	igraph_arpack_options_init(&options);
	info = igraph_authority_score(&graph, &authority_score, NULL, 0, &options);
	if(info!=0)		
	{ puts("ERROR occured during authority score!");
	  exit(-1);
	}

	/* make csv file */
	FILE *csv = fopen(outfile, "w");
	if(graphfile!=NULL)
	{
	  fputs("in.closeness,out.closeness,betweenness,undirected.eigenvector.centrality,pagerank,"
		"hub.score,authority.score\n",
		csv);
	  for(i=0;i<num_vertices;i++)
	  { fprintf(csv, "%g,%g,%g,%g,%g,%g,%g\n", igraph_vector_e(&in_closeness,i),
		    igraph_vector_e(&out_closeness,i), igraph_vector_e(&betweenness,i),
		    absolute(igraph_vector_e(&und_eigen,i)), igraph_vector_e(&pagerank,i),
		    absolute(igraph_vector_e(&hub_score,i)),
		    absolute(igraph_vector_e(&authority_score,i)));
	  }
	}
	fclose(csv);

	/*
		De-allocate objects from memory
	*/
	//igraph_vit_destroy(&vit);  // Vertex iterator
	igraph_vs_destroy(&all_vertices);
	igraph_vector_destroy(&in_closeness);
	igraph_vector_destroy(&out_closeness);
	igraph_vector_destroy(&betweenness);
	igraph_vector_destroy(&und_eigen);
	igraph_vector_destroy(&pagerank);
	igraph_vector_destroy(&hub_score);
	igraph_vector_destroy(&authority_score);
	igraph_destroy(&graph);

	return 0;
}
