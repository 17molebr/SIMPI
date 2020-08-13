#include <unistd.h> 
#include <stdio.h> 
#include <sys/socket.h> 
#include <stdlib.h> 
#include <netinet/in.h> 
#include <string.h> 
#include <iostream>
#include <array>
#include <netdb.h>
#include <thread>
#include <sys/ioctl.h>
#include <signal.h>
#include <arpa/inet.h> 
#include <string.h> 

void new_connection (int sock); 

class server {
    public:
        const char *port= ":8080";
        int make_sock(const char *servspec){
            const int one = 1;
            struct addrinfo hints = {};
            struct addrinfo *res = 0, *ai = 0, *ai4 = 0;
            char *node = strdup(servspec);
            char *service = strrchr(node, ':');
            int sock;

            hints.ai_family = PF_UNSPEC;
            hints.ai_socktype = SOCK_STREAM;
            hints.ai_flags = AI_PASSIVE;

            *service++ = '\0';
            getaddrinfo(*node ? node : "0::0", service, &hints, &res);
            free(node);

            for (ai = res; ai; ai = ai->ai_next) {
                if (ai->ai_family == PF_INET6) break;
                else if (ai->ai_family == PF_INET) ai4 = ai;
            }
            ai = ai ? ai : ai4;

            sock = socket(ai->ai_family, SOCK_STREAM, 0);
            setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
            bind(sock, ai->ai_addr, ai->ai_addrlen);
            listen(sock, 256);
            freeaddrinfo(res);

            return sock;
        }

      
        void accept_loop(const char *servspec){
            int sock = make_sock(servspec);

            for (;;) {
                int new_sock = accept(sock, 0, 0);
                std::thread t(new_connection, new_sock);
                t.detach();
            }
        }

        bool isclosed(int sock) {
            fd_set rfd;
            FD_ZERO(&rfd);
            FD_SET(sock, &rfd);
            timeval tv = { 0 };
            select(sock+1, &rfd, 0, 0, &tv);
            if (!FD_ISSET(sock, &rfd))
                return false;
            int n = 0;
            ioctl(sock, FIONREAD, &n);
            return n == 0;
        }
        
};
void new_connection (int sock) {
     ssize_t r;
     while (!isclosed(sock)) {
        r = send(sock, ".\n", 2, 0);
        if (r < 0) break;
            sleep(1);
        }
        close(sock);
}

int main(int argc, char *argv[]){
    server s;
    signal(SIGPIPE, SIG_IGN);
    s.accept_loop(s.port);
    return 0;
}
