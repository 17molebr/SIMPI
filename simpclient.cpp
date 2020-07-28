#include <stdio.h> 
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <unistd.h> 
#include <string.h> 
#include <stdlib.h>
#include <iostream>
#include <array>
#define PORT 8080 

void run_simpi(int numprocess);  

int main(int argc, char *argv[]) 
{ 
    if (argc != 2) {
        printf("Usage: ./client <num workers>\n");
        exit(2);
    }

    int sock = 0, valread; 
    struct sockaddr_in serv_addr; 
    char *num = argv[1];
    int numprocess = atoi(num);
    if(numprocess <= 4){
        run_simpi(numprocess);
        exit(0); 
    }
    int fullWorkstations;
    int remainderCores; 
    if(numprocess % 4 == 0){
        fullWorkstations = numprocess/4 -1;
        remainderCores = 0;
    }
    else{
        fullWorkstations = numprocess/4;
        remainderCores = numprocess%4;
    }
    char server_ips[8][30] = {"192.168.168.17",
                              "192.168.168.26",
                              "192.168.168.6",
                              "192.168.168.60",
                              "192.168.168.61",
                              "192.168.168.66",
                              "192.168.168.67",
                              "192.168.168.89"}; 
    int numMachines = sizeof(server_ips)/sizeof(server_ips[0]);
    int buffer[1024] = {0}; 
    
    
    //Local Workstaion simpi run
    run_simpi(4); 
    int fullrun = 4;
    for(int i; i < fullWorkstations; i++){   
    
        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
        { 
            printf("\n Socket creation error \n"); 
            return -1; 
        } 
    
        
        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(PORT); 
        // Convert IPv4 and IPv6 addresses from text to binary form 

        if(inet_pton(AF_INET, server_ips[i], &serv_addr.sin_addr)<=0)  
        { 
            printf("\nInvalid address/ Address not supported \n"); 
            return -1; 
        } 
        
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
        { 
            printf("\nConnection Failed \n"); 
            return -1; 
        } 
        if(i == fullWorkstations -1 && remainderCores != 0){
            send(sock , &remainderCores, sizeof(remainderCores), 0 ); 
            send(sock, &numMachines, sizeof(numMachines), 0);
            printf("Nums sent\n");
        }else{
            send(sock , &fullrun, sizeof(fullrun), 0 );
            send(sock, &numMachines, sizeof(numMachines), 0);
            printf("Nums sent\n");
        }

        close(sock);
       
    } 
     
    return 0;
} 

void run_simpi(int numprocess){
    char progname[100];
    strcpy(progname, "mpi");
    char user[100];
    strcpy(user, "user");
    std::string worker_count_str = std::to_string(numprocess);
    char* machines[1];
    strcpy(machines, "4");
    char * args[5] = {progname, user,  const_cast<char*>(worker_count_str.c_str()), machines, NULL};
    if(fork()==0)
        {
            std::cout << progname << "\n";
            std::cout << user << "\n";
            std::cout << numprocess << "\n";
            execv(progname, args);
        }
    return;
}
