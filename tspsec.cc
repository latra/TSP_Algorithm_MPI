/* ******************************************************************** */
/*               Algoritmo Branch-And-Bound Secuencial                  */
/* ******************************************************************** */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <mpi.h>
#include "libtsp.h"
 
using namespace std;
 
unsigned int NCIUDADES;

int main (int argc, char **argv) {
	/*----------------- Inicio de la comunicación MPI -----------------------------*/
	MPI_Status st;
	
	MPI::Init(argc,argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size < 2) {
		printf("Error: Como mínimo deben ejecutarse 2 procesos.\n");
		exit(1);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	switch (argc) {
		case 3:		NCIUDADES = atoi(argv[1]);
					break;
		default:	cerr << "La sintaxis es: bbseq <tamaño> <archivo>" << endl;
					exit(1);
					break;
	}
	/*----------------- Reserva de espacio en memoria e inicialización de variables -----------------------------*/
	int** tsp0 = reservarMatrizCuadrada(NCIUDADES);
	tNodo	nodo,         // nodo a explorar
			lnodo,        // hijo izquierdo
			rnodo;        // hijo derecho
	bool activo,        // condicion de fin
		nueva_U;       // hay nuevo valor de c.s.
	int  U;             // valor de c.s.
	InicNodo (&nodo);              // inicializa estructura nodo

	/*----------------- Lectura de la matriz -----------------------------*/
	if (rank == 0)
		LeerMatriz (argv[2], tsp0);    // lee matriz de fichero
	

	/*----------------- Envío de la matriz a los otros procesos -----------------------------*/
	// Enviamos la array por filas 
	for (int i=0; i<NCIUDADES; i++) {
		MPI_Bcast(tsp0[i], NCIUDADES, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	char buff[10000];
	double t=MPI::Wtime();
	/*----------------- Proceso Mater (rank 0) -----------------------------*/
	int recvType;                  // inicializa cota superior
	if (rank == 0)
	{
		
		/* ------------------- Inicialización de variables exclusivas del proceso Master -------------------*/
		U = INFINITO;
		int slaveStatus[size];
		for (int i = 0; i < size; i++)
			slaveStatus[i] = STATUS_SLEEP;
		tPila pila;         // pila de nodos a explorar
		tNodo solucion;
		PilaInic (&pila);              // inicializa pila
		if (!PilaPush (&pila, &nodo)) {
			printf("Pila llena D:\n");
			liberarMatriz(tsp0);
			exit (1);
		}
		
		while (StillRunning(&pila, slaveStatus, size - 1))
		{

			int position = 0;
			MPI_Recv( buff, 1000, MPI_PACKED, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &st);

			MPI_Unpack(buff, 1000, &position, &recvType, 1, MPI_INT, MPI_COMM_WORLD); 

			switch (recvType)
			{
				// UN proceso eslcavo solicita un trabajo
				case JOB_REQUEST:
					// Revisamos si hay trabajos en cola
					activo = PilaPop (&pila, &nodo);
					if (activo) {
						MPI_Request request;
						position = 0;
						MPI_Pack(&JOB_SENDED, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&nodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&nodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(nodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(nodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						
						// Si había algún trabajo, se lo asignamos al trabajo
						MPI_Isend(buff, 1000, MPI_PACKED, st.MPI_SOURCE,10, MPI_COMM_WORLD, &request);
						slaveStatus[st.MPI_SOURCE - 1] = STATUS_WORKING;
					} else {
						slaveStatus[st.MPI_SOURCE - 1] = STATUS_SLEEP;
					}
					break;
				case NODE_SOLUTION:
				case NODE_PUSH:

					MPI_Unpack(buff, 1000, &position, &nodo.ci, 1, MPI_INT, MPI_COMM_WORLD); 
					if (nodo.ci < U){
						MPI_Unpack(buff, 1000, &position, &nodo.orig_excl, 1, MPI_INT, MPI_COMM_WORLD); 
						MPI_Unpack(buff, 1000, &position, nodo.dest_excl, NCIUDADES -2 , MPI_INT, MPI_COMM_WORLD); 
						MPI_Unpack(buff, 1000, &position, nodo.incl, NCIUDADES, MPI_INT, MPI_COMM_WORLD); 
						if (recvType == NODE_SOLUTION){
							U = nodo.ci;
							CopiaNodo (&nodo, &solucion);
							PilaAcotar (&pila, U);
						}
						else if (recvType == NODE_PUSH) {
							if (!PilaPush (&pila, &nodo)) {
								printf ("Error: pila agotada\n");
								liberarMatriz(tsp0);
								exit (1);
							}
							for (int indice = 0; indice < size - 1 ; indice++) {
								if (slaveStatus[indice] == STATUS_SLEEP) {
									activo = PilaPop (&pila, &nodo);
									if (activo) {
										position = 0;
										MPI_Pack(&JOB_SENDED, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
										MPI_Pack(&nodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
										MPI_Pack(&nodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
										MPI_Pack(nodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
										MPI_Pack(nodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
										
										// Si había algún trabajo, se lo asignamos al trabajo
										MPI_Send(buff, 1000, MPI_PACKED, indice + 1,10, MPI_COMM_WORLD);
										slaveStatus[indice] = STATUS_WORKING;
									}
								}
							}
							// Faltaria ver si hay algun nodo ocioso al que enviarle trabajo
						}
					}
					break;
				default:
					break;
			}
			//EscribeSolucion(&my_message., t);
		}
 		t=MPI::Wtime()-t;
		
		EscribeSolucion(&solucion, t, size);
		liberarMatriz(tsp0);
		printf("Finalizado...\n");
	}
	/* --------------------- Ejecución de los procesos esclavos -----------*/
	else {
		// Solicitud de trabajo
		while(true){
			int position = 0;
			MPI_Pack(&JOB_REQUEST, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
			// Espera de respuesta del master
			MPI_Send( buff, 1000, MPI_PACKED, 0,10, MPI_COMM_WORLD);
			MPI_Recv( buff, 1000, MPI_PACKED, 0, 10, MPI_COMM_WORLD, &st);
			position = 0;
			MPI_Unpack(buff, 1000, &position, &recvType, 1, MPI_INT, MPI_COMM_WORLD); 

			switch (recvType)
			{
				case JOB_SENDED:

					MPI_Unpack(buff, 1000, &position, &nodo.ci, 1, MPI_INT, MPI_COMM_WORLD); 
					MPI_Unpack(buff, 1000, &position, &nodo.orig_excl, 1, MPI_INT, MPI_COMM_WORLD); 
					MPI_Unpack(buff, 1000, &position, nodo.dest_excl, NCIUDADES -2 , MPI_INT, MPI_COMM_WORLD); 
					MPI_Unpack(buff, 1000, &position, nodo.incl, NCIUDADES, MPI_INT, MPI_COMM_WORLD); 

					/* ----------  Inicio tratamiento nodo --------- */
					position = 0;
					Ramifica (&nodo, &lnodo, &rnodo, tsp0);
					int buffer_attach_size = 1000+ MPI_BSEND_OVERHEAD ;
					MPI_Request request;
					if (Solucion(&rnodo)) {
						// Enviar solucion
						position = 0;
						MPI_Pack(&NODE_SOLUTION, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&rnodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&rnodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(rnodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(rnodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						
						MPI_Isend( buff, 1000, MPI_PACKED, 0, 10, MPI_COMM_WORLD, &request);
					}
					else {    
						// Solicitar push
						position = 0;

						MPI_Pack(&NODE_PUSH, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&rnodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&rnodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(rnodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(rnodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Isend( buff, 1000, MPI_PACKED, 0, 10, MPI_COMM_WORLD, &request);
						
					}
					if (Solucion(&lnodo)) {

						position = 0;
						MPI_Pack(&NODE_SOLUTION, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&lnodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&lnodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(lnodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(lnodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						
						MPI_Isend( buff, 1000, MPI_PACKED, 0, 10, MPI_COMM_WORLD,  &request);

					}
					else {   
						// Solicitar push

						position = 0;
						MPI_Pack(&NODE_PUSH, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&lnodo.ci, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(&lnodo.orig_excl, 1, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(lnodo.dest_excl, NCIUDADES - 2, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						MPI_Pack(lnodo.incl, NCIUDADES, MPI_INT, buff, 1000, &position, MPI_COMM_WORLD);
						
						MPI_Isend( buff, 1000, MPI_PACKED, 0, 10, MPI_COMM_WORLD,  &request);


					}
					break;
			}
		}
	}

}