#include "simpi.h"
//FIX server hard code in new_connection around 254
//Fix hardcode in new_connection status 3 case to num workstations 
// init global simpi
static simpi *main_simpi;
struct timeval begin, end1;

// simpi init function
void SIMPI_INIT(int par_id, size_t synch_size, int workstations, int workstationid)
{
    main_simpi = new simpi(par_id, synch_size, workstations, workstationid);
}
// simpi synch
void SIMPI_SYNCH()
{
    main_simpi->synch();
}

void SIMPI_FINALIZE()
{
    //main_simpi->synch();
    delete main_simpi;
    char buff[100] = "pkill -f user";
    system(buff);
    return;
}

class server;
void new_connection(int sock, server s);
void new_connection2(int sock, server s);
std::vector<int> workstations;



void run_client(matrix m, int s){
    /*
    data_info info;
    info.arr = m.arr;
    */
    //std::cout << s << "\n";
    /*
    char array[2]; 
    int r = read(s, &array, 2);
    if(r < 0){
        printf("read error\n");
    }
    std::cout << array; 
    */
    ssize_t r;
    int sendval = 1;
    int start = main_simpi->get_start();
    int end = main_simpi->get_end();
    int id = main_simpi->get_workstation_id();
    int xdim = m.get_x();
    int ydim = m.get_y();
    r = send(s, &sendval, sizeof(sendval), 0);
    r = send(s, &start, sizeof(start), 0); 
    r = send(s, &end, sizeof(end), 0); 
    r = send(s, &xdim, sizeof(xdim), 0); 
    r = send(s, &ydim, sizeof(ydim), 0);
    r = send(s, &id, sizeof(id), 0);
    std::cout << m.arr[6*xdim + 4];
    for (int a = start; a < end; a++)
    {
        for (int b = 0; b < xdim; b++)
        {
            double element = m.arr[a*xdim + b];
            //std::cout << element << "\n";
            r = send(s, &element, sizeof(element), 0);
        }
    }
    //close(s);
    return;
}

int run_client2(matrix &m, int s){
    int r;
    int sendval = 2;
    int done;
    int id = main_simpi->get_workstation_id();
    int xdim = m.get_x();
    int ydim = m.get_y();
    r = send(s, &sendval, sizeof(sendval), 0);
    r = send(s, &id, sizeof(id), 0);
    r = read(s, &done, sizeof(done));
    if(done == 0){
        return 0;
    }
    std::cout << "\nin distribute done = \n"<< done;
    for (int a = 0; a < ydim; a++){
        for (int b = 0; b < xdim; b++){
            double element = 0;
            r = read(s, &element, sizeof(element));
            m.arr[a*10 +b] = element;
            //std::cout << "\nElement in client is"<< element << "\n";
        }
    }
    //close(s);
    return 1; 
}

void run_client_synch(int s){
    //synch function to wait for all parts of matix to be completed
    int status = 3;
    int client_status = 0;
    int value = 0;
    client_status = send(s, &status, sizeof(status), 0);
    while(value == 0){
        client_status = read(s, &value, sizeof(value));
        std::cout << "\nstaus = " << client_status <<"\n";
    }
    return;
}

class server {
    public:
        int count = 0; //Number of machines has accessed 
        const char *port= ":8080";
        int defualt_size = 50;
        int num_runs = 0;
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
            if(getaddrinfo(*node ? node : "0::0", service, &hints, &res)!=0){
                printf("Error with getaddrinfo\n");
            }
            free(node);

            for (ai = res; ai; ai = ai->ai_next) {
                if (ai->ai_family == PF_INET6) break;
                else if (ai->ai_family == PF_INET) ai4 = ai;
            }
            ai = ai ? ai : ai4;

            sock = socket(ai->ai_family, SOCK_STREAM, 0);
            if(sock < 0){
                printf("Socket creation erorr\n");
            }
            if(setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one))!= 0){
                printf("sockoptions erorr");
            }
            if(bind(sock, ai->ai_addr, ai->ai_addrlen)!=0){
                perror("bind erorr: ");
            }
            printf("sock = : %d",sock);
            listen(sock, 256);
            freeaddrinfo(res);

            return sock;
        }

        void accept_loop(const char *servspec, server s){
            int sock = make_sock(servspec);
            
            for (;;) {
                int new_sock = accept(sock, 0, 0);
                std::thread t(new_connection, new_sock, s);
                t.detach();
                /*
                count += 1;
                if(count == 1){
                    count = 0;
                    return;
                }
                */
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

client c;
server s;
double *temp= new double[s.defualt_size];
int num = 0;
int *num_connections = &num;
int *workstation_status = new int[2];
int current_x = 0;
int current_y = 0;

void new_connection(int sock, server s) {
    //Server Connection 
    ssize_t r;
    /*
    data_info info;
    */
    while (!s.isclosed(sock)) {
        /*
        r = send(sock, ".\n", 2, 0);
        if (r < 0) 
            break;
        sleep(1)*/
        /*
        int workstaion_id = 0;
        int stage = 0;
        r = read(sock, &workstaion_id, sizeof(workstaion_id));
        stage = workstations.at(workstaion_id-1);
        */

        

        //
        int status = 0;
        r = read(sock, &status, sizeof(status));
        std::cout << "\nStatus = " << status << "\n"; 
        if(status == 1){
            int start;
            int end;
            int xdim;
            int ydim;
            int id;
            r = read(sock, &start, sizeof(start));
            r = read(sock, &end, sizeof(end));
            std::cout << start << "\n";
            std::cout << end << "\n";
            r = read(sock, &xdim, sizeof(xdim));
            r = read(sock, &ydim, sizeof(ydim));
            r = read(sock, &id, sizeof(id));
            if(xdim * ydim != s.defualt_size && s.num_runs == 0 && start == 0){
                double *array = new double[xdim*ydim];
                delete [] temp;
                temp = array;
                current_x = xdim;
                current_y = ydim;
            }
            //must synch after this step to make sure that you are not attempting to write to something that has been deleated 
            //*synch*
            workstation_status[id] = 1;
            while(1){
                int flag = 0;
                for(int i =1; i <= 2; i++){ //Fix this to num_workstations 
                    if(workstation_status[i] == 0){
                        flag = 1;
                    }
                }
                if(flag != 1){
                    break;
                }
            }
            workstation_status[id] = 0;
            //*end synch*
            for (int a = start; a < end; a++)
            {
                for (int b = 0; b < xdim; b++)
                {
                    double element = 0;
                    r = read(sock, &element, sizeof(element));
                    temp[a*xdim +b] = element;
                    //std::cout << elemeent << "\n";
                }
            }
    
            
            for (int i = 0; i < xdim; i++)
            {
                std::cout << "\n";
                for (int j = 0; j < ydim; j++)
                {
                    std::cout << std::fixed << std::setprecision(2) << temp[i + j * xdim];
                    std::cout << ", ";
                }
            }
            std::cout << "\n";
            s.num_runs += 1;
            std::cout <<"\nIncrementing ip = " << id <<"\n";
            workstation_status[id] = 1;

        }
        if(status == 2){
            //check for matrix completion if not then send back to client that unable to disribute yet 
            std::cout << "\n" << "Server in status = 2" << "\n";
            std::cout << "\nConnection 1 = " << workstation_status[1] << "\n";
            std::cout << "\nConnection 2 = " << workstation_status[2] << "\n";
            int status1 = 0;
            int send2 = 1;
            int id = 0;
            status1 = read(sock, &id, sizeof(id));
            status1 = send(sock, &send2, sizeof(send2), 0); //ok to distribute 
            std::cout << "\n" << "outside if statement" << "\n";
            s.num_runs = 0;
            std::cout << "\n" << "Matrix has been redistributed"<<"\n";
            std::cout << "\n x = " << current_x << " y = " << current_y << "\n";
            for (int a = 0; a < current_y; a++){
                for (int b = 0; b < current_x; b++){
                    double element = temp[a*current_x + b];
                    //std::cout << "\nElement is "<< element << "\n";
                    r = send(sock, &element, sizeof(element), 0);
                }
            }
            std::cout << "\n" << "Matrix has been redistributed"<<"\n";
            workstation_status[id] = 0;
            //clear workstation_status value if last distribute
        }
        if(status == 3){
            //Wait for whole matix to be updated and then send to client to move on and then start redistributing
            while(1){
                int flag = 0;
                for(int i =1; i <= 2; i++){ //Fix this to num_workstations 
                    if(workstation_status[i] == 0){
                        flag = 1;
                    }
                }
                if(flag != 1){
                    break;
                }
            }
            int status2 = 0;
            int send3 = 1; //Send value to have client break wating loop
            status2 = send(sock, &send3, sizeof(send3), 0);
        }
    }
    
    close(sock);
    
}


/******************Simpi Functions*************************/
simpi::simpi(int _id, int _num_workers, int _num_workstaions, int _workstation_id)
{
    id = _id;
    num_workers = _num_workers;
    num_workstations = _num_workstaions;
    workstationid = _workstation_id;
    size_t synchObjectSize =
        sizeof(synch_object) + sizeof(int) * (num_workers + 1);
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
    if (fd == -1)
    {
        perror("Unable to shm_open synch_object: ");
        exit(1);
    }
    shm_fd = fd;
    synch_info = (synch_object *)mmap(NULL, synchObjectSize,
                                      PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (synch_info == MAP_FAILED)
    {
        perror("Unable to mmap synch_info: ");
        exit(1);
    }
    std::cout << "workstation id = " << workstationid <<" " << id;
    if(workstationid == 0 && id == 0){
        std::cout << "serversetup";
        gettimeofday(&begin, 0);
        workstations.reserve(num_workstations);
        std::fill(workstations.begin(), workstations.end(), 0);
        //signal(SIGPIPE, SIG_IGN);
        std::cout << "signal sent, launching accept loop";
        s.accept_loop(s.port, s);
        //s.sendback_accept_loop(s.port, s);
        return;
    }
    else if(workstationid == 0){
        std::cout << "workstationid == 0, exiting";
        exit(0);
    }
    else if(workstationid != 0 && id == 0){
        std::cout << "client setup workstation id =" << workstationid;
        c.setup_client();
        //std::cout << c.sock;
    }
       
}
simpi::~simpi()
{
    for (std::pair<std::string, matrix_metadata> matrix : matrix_info)
    {
        free_matrix(matrix.second.unique_id);
    }
    shm_unlink(SYNCH_OBJECT_MEM_NAME);
    close(shm_fd);
}
void simpi::synch()
{
    int *ready = synch_info->ready;
    int synchid = ready[num_workers] + 1;
    ready[id] = synchid;
    int breakout = 0;
    while (1)
    {
        breakout = 1;
        for (int i = 0; i < num_workers; i++)
        {
            if (ready[i] < synchid)
            {
                breakout = 0;
                break;
            }
        }
        if (breakout == 1)
        {
            ready[num_workers] = synchid;
            // and here we increment the additional variable
            break;
        }
    }
}

void SIMPI_DISTRIBUTE(matrix m, matrix &m1){
    if(main_simpi->get_workstation_id() != 0 && main_simpi->get_id()== 0){
        std::cout << "passed in sock= "<< c.sock << "\n";
        run_client(m, c.sock);
        int done_val;
        run_client_synch(c.sock);
        done_val = run_client2(m1, c.sock);
        return;
    }
    return;
}

std::pair<std::string, double *> simpi::create_matrix(int x, int y)
{
    size_t size = x * y * sizeof(double);

    if (id == 0)
    {
        // generate a uniqueid for matrix
        std::string unique_id = get_shared_mem_name();
        // create a shared mem object
        int fd = shm_open(unique_id.c_str(), O_RDWR | O_CREAT, 0777);
        if (fd == -1)
        {
            std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                              unique_id + "):  ";
            perror(msg.c_str());
            exit(1);
        }
        ftruncate(fd, size);

        // allocate matrix
        double *matrix =
            (double *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (matrix == MAP_FAILED)
        {
            perror("Unable to mmap matrix: ");
            exit(1);
        }

        matrix_metadata metadata;
        metadata.size = size;
        metadata.file_descriptor = fd;
        strcpy(metadata.unique_id, unique_id.c_str());
        metadata.matrix_data = matrix;
        matrix_info[unique_id] = metadata;

        // write name to synch_object so that other processes can get the uniqie
        // id
        strcpy(synch_info->last_matrix_id, unique_id.c_str());
        synch();
        std::pair<std::string, double *> pass_back;
        pass_back = std::make_pair(unique_id, matrix);
        return pass_back;
    }
    else
    {
        // wait for id 0 to create the shared memory for matrix
        synch();
        // get the unique id from the synch object
        std::string unique_id = synch_info->last_matrix_id;
        // open and allocate the shared memory
        int fd = shm_open(unique_id.c_str(), O_RDWR, 0777);
        if (fd == -1)
        {
            std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                              unique_id + "):  ";
            perror(msg.c_str());
            exit(1);
        }
        double *matrix =
            (double *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (matrix == MAP_FAILED)
        {
            std::string msg = std::to_string(id) + ": Unable to mmap matrix: ";
            perror(msg.c_str());
            exit(1);
        }

        // create a metadata object
        matrix_metadata metadata;
        strcpy(metadata.unique_id, unique_id.c_str());
        metadata.file_descriptor = fd;
        metadata.matrix_data = matrix;
        matrix_info[unique_id] = metadata;
        std::pair<std::string, double *> pass_back =
            std::make_pair(unique_id, matrix);
        return pass_back;
    }
}

void simpi::free_matrix(std::string unique_id)
{
    // get matrix metadata
    // close fd, shmunmap, munmap
    matrix_metadata metadata = matrix_info[unique_id];
    close(metadata.file_descriptor);
    shm_unlink(metadata.unique_id);
    munmap(metadata.matrix_data, metadata.size);
}

std::string simpi::get_shared_mem_name()
{
    // gets a unique name for each shared memory based on the time that each was
    // made
    timeval curTime;
    gettimeofday(&curTime, NULL);
    unsigned long micro = curTime.tv_sec * (uint64_t)1000000 + curTime.tv_usec;
    std::string name = "simpi_" + std::to_string(micro);
    return name;
}


/******************Matrix Functions*************************/
matrix::matrix(int x, int y) // constructor
{
    // use main_simp init the matrix for all processes. The id is also in simp
    std::pair<std::string, double *> pass_back(main_simpi->create_matrix(x, y));
    unique_id = pass_back.first;
    arr = pass_back.second;
    xdim = x;
    ydim = y;
}
matrix::~matrix() // destructor
{
    // use main_simpi for getting rid of the mem and unlink stuff
    main_simpi->free_matrix(unique_id);
}

std::ostream &operator<<(std::ostream &out, const matrix &m)
{
    if (main_simpi->get_id() == 0)
    {
        for (int i = 0; i < m.xdim; i++)
        {
            out << "\n";
            for (int j = 0; j < m.ydim; j++)
            {
                out << std::fixed << std::setprecision(2) << m.arr[i + j * m.xdim];
                out << ", ";
            }
        }
        out << "\n";
        return out;
    }
    else
    {
        return out;
    }
}
/*
Old and slower version of finding matrix inverse
Stored as an archive
*/
void matrix::inverse_old(matrix *inverse)
{
    if (get_x() != get_y())
    {
        std::cout << "Invalid Matrix";
        exit(1);
    }

    matrix *adj = new matrix(get_x(), get_y());

    // Find determinant of A[][]
    int det = determinant(arr, get_x(), get_y());
    if (det == 0)
    {
        std::cout << "Singular matrix, can't find its inverse";
        return;
    }
    main_simpi->synch();

    adjoint(arr, adj->arr, xdim, main_simpi->get_id(),
            main_simpi->get_num_workers());
    main_simpi->synch();

    int num_processes = main_simpi->get_num_workers();
    int parID = main_simpi->get_id();
    int n = get_x();
    if (num_processes > n)
    {
        num_processes = n;
    }

    int rpp = n / num_processes;
    int start = main_simpi->get_id() * rpp;
    int end = start + rpp;

    for (int i = start; i < end; i++)
        for (int j = 0; j < get_x(); j++)
        {
            inverse->get(i, j) = (adj->get(i, j)) / double(det);
        }
    if (n % num_processes != 0)
    {
        int leftover = n % num_processes;
        if (parID < leftover)
        {
            parID += (n - leftover);
            int start = parID;
            int end = start + 1;
            for (int i = start; i < end; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverse->get(i, j) = (adj->get(i, j)) / double(det);
                }
            }
        }
    }
    main_simpi->synch();
    return;
}

/*
This is a helper function to calculate the determinant of a matrix
*/
int matrix::determinant(double *A, int n, int order)
{
    int D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0];

    double temp[order * order]; // To store cofactors

    int sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n, order);
        D += sign * A[0 + f * order] * determinant(temp, n - 1, order);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

/*
This is a helper function to calculate the adjoint of a matrix
*/
void matrix::adjoint(double *A,
                     double *adj,
                     int order,
                     int par_id,
                     int par_count)
{
    if (order == 1)
    {
        adj[0] = 1;
        return;
    }

    int rpp = order / par_count;
    int start = par_id * rpp;
    int end = start + rpp;

    // temp is used to store cofactors of A[][]
    int sign = 1;
    double temp[order * order];

    for (int i = 0; i < order; i++)
    {
        for (int j = start; j < end; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, order, order);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j + i * order] = (sign) * (determinant(temp, order - 1, order));
        }
    }
}

/*
This is a helper function to calculate the cofactor of a matrix
*/
void matrix::getCofactor(double *A,
                         double *temp,
                         int p,
                         int q,
                         int n,
                         int order)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[(i) + (j++) * order] = A[row + col * order];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/*
This function calculates the lower and upper triangular matrices of a square nxn matrix
It requires an input of 2 empty nxn matrices that are modified
A = LU
*/
void matrix::luDecomposition(matrix *lower, matrix *upper)
{

    // Check if Matrix is square
    if (get_x() != get_y())
    {
        std::cout << "Invalid Matrix";
        exit(1);
    }

    for (int i = 0; i < get_x(); i++)
    {
        // Calculate work per parallel process
        // Has to be calculated on every loop iteration as the inner loop is decrementing
        int num_processes = main_simpi->get_num_workers();
        int parID = main_simpi->get_id();
        int total = xdim - i;
        if (num_processes > total)
        {
            num_processes = total;
        }
        int rpp = total / num_processes;
        int start = rpp * main_simpi->get_id() + i;
        int end = start + rpp;

        // Upper Triangular
        for (int k = start; k < end; k++)
        {
            if (k >= get_x())
            {
                break;
            }
            // Summation of L(i, j) * U(j, k)
            float sum = 0;
            for (int j = 0; j < i; j++)
                sum += (lower->get(i, j) * upper->get(j, k));

            // Evaluating U(i, k)
            upper->get(i, k) = get(i, k) - sum;
        }

        // Calculate and execute which processes take the leftover work
        if (total % num_processes != 0)
        {
            int leftover = total % num_processes;
            if (parID < leftover)
            {
                parID += (xdim - leftover);
                int start = parID;
                int end = start + 1;
                for (int a = start; a < end; a++)
                {
                    // Summation of L(i, j) * U(j, k)
                    float sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (lower->get(i, j) * upper->get(j, a));
                    // Evaluating U(i, k)
                    upper->get(i, a) = get(i, a) - sum;
                }
            }
        }

        main_simpi->synch();

        total = get_x() - i;
        parID = main_simpi->get_id();
        num_processes = main_simpi->get_num_workers();

        // Lower Triangular
        for (int k = start; k < end; k++)
        {
            if (k >= get_x())
            {
                break;
            }
            if (i == k)
            {
                lower->get(i, i) = 1; // Diagonal as 1
            }
            else
            {
                // Summation of L(k, j) * U(j, i)
                float sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (lower->get(k, j) * upper->get(j, i));
                // Evaluating L(k, i)
                lower->get(k, i) = ((get(k, i) - sum) / upper->get(i, i));
            }
        }

        // Calculate and execute which processes take the leftover work
        if (total % num_processes != 0)
        {
            int leftover = total % num_processes;
            if (parID < leftover)
            {
                parID += (get_x() - leftover);
                int start = parID;
                int end = start + 1;
                for (int a = start; a < end; a++)
                {
                    if (i == a)
                        lower->get(i, i) = 1; // Diagonal as 1
                    else
                    {
                        // Summation of L(k, j) * U(j, i)
                        float sum = 0;
                        for (int j = 0; j < i; j++)
                            sum += (lower->get(a, j) * upper->get(j, i));

                        // Evaluating L(k, i)
                        lower->get(a, i) = (get(a, i) - sum) / upper->get(i, i);
                    }
                }
            }
        }
        main_simpi->synch();
    }
    return;
}

/*
This method calculates the inverse of a matrix by using its LU Decomposition
A = LU
LZ = B
LX = Z
B corresponds to individual columns of an nxn identity matrix
X represents each corresponding column of the inverse matrix
*/
void matrix::inverse(matrix *inv)
{

    //Check if matrix is square
    if (get_x() != get_y())
    {
        std::cout << "Invalid Matrix";
        exit(1);
    }

    //Solve for lower and upper matrices
    matrix *upper = new matrix(get_x(), get_y());
    matrix *lower = new matrix(get_x(), get_y());
    luDecomposition(lower, upper);
    main_simpi->synch();

    //Create Identity nxn Matrix
    matrix *identity = new matrix(get_x(), get_y());
    if (main_simpi->get_id() == 0)
    {
        for (int i = 0; i < get_x(); i++)
        {
            for (int j = 0; j < get_x(); j++)
            {
                if (i == j)
                {
                    identity->get(i, j) = 1;
                }
                else
                {
                    identity->get(i, j) = 0;
                }
            }
        }
    }

    main_simpi->synch();

    // Calculate columns per parallel process
    int num_workers = main_simpi->get_num_workers();
    if (num_workers > get_x())
    {
        num_workers = get_x();
    }
    int cpp = get_x() / num_workers;
    int start = cpp * main_simpi->get_id();
    int end = start + cpp;

    // Initialize necessary arrays for future calculations
    // Each array is local to its own process
    float identity_col[get_x()];
    float z_col[get_x()];
    float soln_col[get_x()];

    for (int a = start; a < end; a++)
    {

        //Get individual columns of identity matrix
        for (int b = 0; b < get_x(); b++)
        {
            identity_col[b] = identity->get(b, a);
        }

        //Reset Z column to solve for again
        for (int d = 0; d < get_x(); d++)
        {
            z_col[d] = 0;
        }

        //Solve LZ = I
        (*lower).forward_substitution(identity_col, z_col);

        //Reset X column to solve for again
        for (int d = 0; d < get_x(); d++)
        {
            soln_col[d] = 0;
        }

        //Solve UX = Z
        (*upper).backward_substitution(z_col, soln_col);

        //Input X column to corresponding columnn in final inverse matrix
        for (int c = 0; c < get_x(); c++)
        {
            inv->get(c, a) = soln_col[c];
        }
    }

    // Calculate and execute which processes take the leftover rows
    // ex. with 3 processes and a 10x10 matrix:
    // 0-3 is process 0
    // 3-6 is process 1
    // 6-9 is process 2
    // 9 is the leftover column that is taken by process 0
    int parID = main_simpi->get_id();
    if (get_x() % num_workers != 0)
    {
        int leftover = get_x() % num_workers;
        if (parID < leftover)
        {
            parID += (get_x() - leftover);
            int start = parID;
            int end = start + 1;
            for (int a = start; a < end; a++)
            {
                //Get individual columns of identity matrix
                for (int b = 0; b < get_x(); b++)
                {
                    identity_col[b] = identity->get(b, a);
                }

                //Reset Z column to solve for again
                for (int d = 0; d < get_x(); d++)
                {
                    z_col[d] = 0;
                }

                //Solve LZ = I
                (*lower).forward_substitution(identity_col, z_col);

                //Reset X column to solve for again
                for (int d = 0; d < get_x(); d++)
                {
                    soln_col[d] = 0;
                }

                //Solve UX = Z
                (*upper).backward_substitution(z_col, soln_col);

                //Input X column to corresponding columnn in final inverse matrix
                for (int c = 0; c < get_x(); c++)
                {
                    inv->get(c, a) = soln_col[c];
                }
            }
        }
    }

    main_simpi->synch();
    return;
}

/*
This is a helper function to calculate the solutions of a lower triangular matrix
*/
void matrix::forward_substitution(float *b, float *x)
{
    double suma;
    for (int i = 0; i < get_x(); i = i + 1)
    {
        suma = 0;
        for (int j = 0; j < i; j = j + 1)
            suma = suma + get(i, j) * x[j];

        x[i] = (b[i] - suma) / get(i, i);
    }
}

/*
This is a helper function to calculate the solutions of an upper triangular matrix
*/
void matrix::backward_substitution(float *b, float *x)
{
    double suma;
    for (int i = get_x() - 1; i >= 0; i = i - 1)
    {
        suma = 0;
        for (int j = get_x() - 1; j > i; j = j - 1)
            suma = suma + get(i, j) * x[j];

        x[i] = (b[i] - suma) / get(i, i);
    }
}

/*
Solves a linear system of equations in parallel if the matrix is diagonally dominant
outputs solution to a vector of 0s passed into it
void->void
*/

void matrix::jacobi(vector *constants, vector *solution)
{

    int processCount = main_simpi->get_num_workers();
    int id = main_simpi->get_id();

    vector *prev = new vector(constants->get_size()); // shared mem containing a copy of values
    //vector *solution = new vector(*main_simpi, constants->get_size()); // shared mem containing actual calculated values

    matrix *saveEq = new matrix(get_x(), get_y() + 1);     // save equations from modification
    vector *saveConst = new vector(constants->get_size()); // saves input vector

    int work = constants->get_size() / processCount;

    int i, j, k;
    int start = id * work;
    int end = start + work;

    main_simpi->synch();

    //Save Matrix and Vector
    for (i = start; i < end; i++)
    {
        for (j = 0; j < get_y(); j++)
        {
            saveEq->get(i, j) = get(i, j);
        }
        saveEq->get(i, get_y() + 1) = constants->get(i);
        saveConst->get(i) = constants->get(i);
    }

    //synch, wait for all process before solving
    main_simpi->synch();

    //setup, switch var coefficient with row solution, and divide by coefficient
    for (i = start; i < end; i++)
    {
        double temp = get(i, i);
        get(i, i) = constants->get(i);
        constants->get(i) = temp;
        for (j = 0; j < get_y(); j++)
        {
            if (j != i)
            {
                get(i, j) *= -1;
            }
            get(i, j) /= constants->get(i);
            //main_simpi->synch();
        }
        prev->get(i) = 1.0;
        solution->get(i) = constants->get(i);
        main_simpi->synch();
    }

    main_simpi->synch();
    // first iteration by trying substituting 1
    for (i = start; i < end; i++)
    {
        double rowSum = 0;
        for (j = 0; j < get_y(); j++)
        {
            if (j == i)
            {
                rowSum += get(i, j);
            }
            else
            {
                rowSum += (get(i, j) * prev->get(j));
            }
            //main_simpi->synch();
        }
        solution->get(i) = rowSum;
        main_simpi->synch();
    }

    //wait for all processes before repeating iterations with calculated results
    //main_simpi->synch();

    for (k = 0; k < 1000; k++)
    {
        for (i = start; i < end; i++)
        {
            //save prev value for comparision
            prev->get(i) = solution->get(i);
            main_simpi->synch();
        }
        for (i = start; i < end; i++)
        {
            double rowSum = 0;
            for (j = 0; j < get_y(); j++)
            {
                if (j == i)
                {
                    rowSum += get(i, j);
                }
                else
                {
                    rowSum += (get(i, j) * prev->get(j));
                }
            }
            solution->get(i) = rowSum;
            //main_simpi->synch();
        }
        //wait at end of each loop for all processes before beginning next iteration
        main_simpi->synch();
    }
    main_simpi->synch();
    //restore original matrix and vector
    for (i = start; i < end; i++)
    {
        for (j = 0; j < get_y(); j++)
        {
            get(i, j) = saveEq->get(i, j);
        }
        constants->get(i) = saveConst->get(i);
    }
    //wait for all processes before returning solution vector
    main_simpi->synch();
    //return solution;
    return;
}

/*
Checks if a square matrix is diagonally dominant
(diagonal terms are greater-than or equal to sum of the rest of their row)
none->bool
*/
bool matrix::isDiagonallyDominant()
{
    for (int i = 0; i < get_x(); i++)
    {
        double sq;
        double rest = 0;
        for (int j = 0; j < get_y(); j++)
        {
            if (i == j)
            {
                sq = get(i, j);
            }
            else
            {
                rest += get(i, j);
            }
        }
        if (sq < rest)
        {
            return false;
        }
    }
    return true;
}

/*
Solves a linear system of equations
If the matrix is diagonally dominant, jacobi-iterative method is used
else, the inverse mutliplication method is used
takes in a vector of constants that each row is to be solved for and a solution vector of 0s in which the solution will be written in
void -> void
*/
void matrix::solveSystem(vector *constants, vector *solution)
{
    bool dd = isDiagonallyDominant();
    main_simpi->synch();
    if (dd)
    {
        main_simpi->synch();
        jacobi(constants, solution);
    }
    else
    {
        if (main_simpi->get_id() == 0)
        {
            std::cout << "failsafe" << std::endl;
        }
        main_simpi->synch();
        failSafe(constants, solution);
    }
    main_simpi->synch();
}

/*
Method to solve a system of linear equations if the system is not diagonally dominant
Uses the inverse-mutliplication method.
void->void
*/
void matrix::failSafe(vector *constants, vector *solution)
{
    matrix *inv = new matrix(get_x(), get_y());
    inverse(inv);
    std::cout << "inverse calculated" << std::endl;
    main_simpi->synch();

    int processCount = main_simpi->get_num_workers();
    int id = main_simpi->get_id();
    int work = constants->get_size() / processCount;

    int start = id * work;
    int end = start + work;
    int n = constants->get_size();

    double sol;
    for (int i = start; i < end; i++)
    {
        sol = 0;
        for (int j = 0; j < n; j++)
        {
            sol += (inv->get(i, j) * constants->get(j));
        }
        solution->get(i) = sol;
    }
    main_simpi->synch();
    return;
}

/******************Vector Functions*************************/
vector::vector(int a)
{
    // use simp and init the matrix for all processes. The id is also in simp
    std::pair<std::string, double *> pass_back(main_simpi->create_matrix(1, a));
    unique_id = pass_back.first;
    arr = pass_back.second;
    dim = a;
}
vector::~vector() // destructor
{
    // use main_simpi for getting rid of the mem and unlink stuff
    main_simpi->free_matrix(unique_id);
}

//ALGEBRA FUNCTIONS :-

/*
Solves matrix multiplication in parallel and outputs the product solution
matrix -> matrix
*/
matrix &matrix::multiply(matrix other)
{
    matrix *result = new matrix(xdim, other.get_y());
    int number_of_processes = main_simpi->get_synch_info()->par_count;
    int number_of_workstations = main_simpi->get_num_workstations();
    std::cout << "Number of processes = " << number_of_processes << std::endl;
    std::cout << "Number of workstations = " << number_of_workstations << std::endl;
    int tempForProcesses = number_of_processes;
    number_of_processes = number_of_processes * number_of_workstations;
    int parId = main_simpi->get_id();
    int workstationid = main_simpi->get_workstation_id() - 1 ;
    int parIDInit = parId;
    parId = workstationid * 4 + parId;

    printf("WORKSTATION ID = %d\n", workstationid);
    printf("parID ID = %d\n", parId);
    int numCols = other.get_y() / number_of_workstations;
    if (parId <= other.get_y())
    {
        if (number_of_processes > other.get_y())
        {
            number_of_processes = other.get_y();
        }
        int Arow = get_x();
        int Acol = get_y();
        int Brow = other.get_x();
        int Bcol = other.get_y();


        if (Acol != Brow)
        {
            printf("Error. Matrices can't be multiplied\n");
        }

        
        int rpp = Bcol / 10;
        int start = rpp * parIDInit + numCols * workstationid;
        int end = start + rpp;
        main_simpi->set_start(workstationid * numCols);
        main_simpi->set_end((workstationid + 1) * numCols);
        if (numCols % tempForProcesses != 0){
            int leftover = numCols % tempForProcesses;
            printf("DEBUG 1\n");
            if (parIDInit < leftover)
                {
                    printf("DEBUG 2 par ID : %d numcols : %d\n", parId, numCols);
                // parId += (Arow - leftover);
                int start =  numCols * workstationid + (numCols - 1 );
                int end = start + 1;
                // main_simpi->set_start(start);
                // main_simpi->set_end(start + (end-start) * tempForProcesses);
                for (int a = start; a < end; a++)
                {
                    for (int b = 0; b < Arow; b++)
                    {
                        int sum = 0;
                        for (int c = 0; c < Brow; c++)
                        {
                            sum = sum + other.get_algbera(c + a * Brow) * get_algbera(c * Arow + b);
                        }
                        result->set(a * Arow + b, sum);
                    }
                }
            }
        }
        for (int a = start; a < end; a++)
        {
            for (int b = 0; b < Arow; b++)
            {
                int sum = 0;
                for (int c = 0; c < Brow; c++)
                {
                    sum = sum + other.get_algbera(c + a * Brow) * get_algbera(c * Arow + b);
                }
                result->set(a * Arow + b, sum);
            }
        }
        main_simpi->synch();
        return *result;
    }
}

/*
Outputs the trasnpose solution of a matrix entered
void -> matrix
*/
matrix &matrix::transpose()
{
    matrix *result = new matrix(get_y(), get_x());
    int number_of_processes = main_simpi->get_synch_info()->par_count;

    if (number_of_processes > get_x())
    {
        number_of_processes = get_x();
    }

    int Arow = get_x();
    int Acol = get_y();

    int rpp = Arow / number_of_processes;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Acol; j++)
        {
            result->set(j + i * Acol, get_algbera(j * Arow + i));
        }
    }
    return *result;
}

/*
Solves matrix multiplication in parallel and outputs the product solution
matrix -> matrix
*/
matrix &matrix::add(matrix other)
{
    matrix *result = new matrix(xdim, other.get_y());
    int number_of_processes = main_simpi->get_synch_info()->par_count;
    int number_of_workstations = main_simpi->get_num_workstations();
    std::cout << "Number of processes = " << number_of_processes << std::endl;
    std::cout << "Number of workstations = " << number_of_workstations << std::endl;
    int tempForProcesses = number_of_processes;
    number_of_processes = number_of_processes * number_of_workstations;
    int parId = main_simpi->get_id();
    int workstationid = main_simpi->get_workstation_id() - 1 ;
    int parIDInit = parId;
    parId = workstationid * 4 + parId;

    printf("WORKSTATION ID = %d\n", workstationid);
    printf("parID ID = %d\n", parId);
    int numCols = other.get_y() / number_of_workstations;
    if (parId <= other.get_y())
    {
        if (number_of_processes > other.get_y())
        {
            number_of_processes = other.get_y();
        }
        int Arow = get_x();
        int Acol = get_y();
        int Brow = other.get_x();
        int Bcol = other.get_y();


        if (Acol != Brow)
        {
            printf("Error. Matrices can't be multiplied\n");
        }

        
        int rpp = Bcol / 10;
        int start = rpp * parIDInit + numCols * workstationid;
        int end = start + rpp;
        main_simpi->set_start(workstationid * numCols);
        main_simpi->set_end((workstationid + 1) * numCols);
        if (numCols % tempForProcesses != 0){
            int leftover = numCols % tempForProcesses;
            printf("DEBUG 1\n");
            if (parIDInit < leftover)
                {
                    printf("DEBUG 2 par ID : %d numcols : %d\n", parId, numCols);
                int start =  numCols * workstationid + (numCols - 1 );
                int end = start + 1;
                for (int a = start; a < end; a++)
                {
                    for (int b = 0; b < Arow; b++)
                    {
                        result->set(b + a * Arow,
                         (other.get_algbera(b + a * Arow) + get_algbera(b + a * Arow)));
                    }
                }
            }
        }
        for (int a = start; a < end; a++)
        {
            for (int b = 0; b < Arow; b++)
            {
                result->set(b + a * Arow,
                         (other.get_algbera(b + a * Arow) + get_algbera(b + a * Arow)));
            }
        }
        main_simpi->synch();
        return *result;
    }
}







// matrix &matrix::add(matrix other)
// {
//     matrix *result = new matrix(xdim, other.get_y());
//     int number_of_processes = main_simpi->get_synch_info()->par_count;
//     int number_of_workstations = main_simpi->get_num_workstations();
//     int tempForProcesses = number_of_processes;
//     number_of_processes = number_of_processes * number_of_workstations;
//     int parId = main_simpi->get_id();
//     int workstationid = main_simpi->get_workstation_id() - 1 ;
//     parId = workstationid * tempForProcesses + parId;
//     printf("WORKSTATION ID = %d\n", workstationid);
//     printf("parID ID = %d\n", parId);
//     if (parId <= other.get_y())
//     {
//         if (number_of_processes > other.get_y())
//         {
//             number_of_processes = other.get_y();
//         }
//         int Arow = get_x();
//         int Acol = get_y();
//         int Brow = other.get_x();
//         int Bcol = other.get_y();
//         int rpp = Bcol / number_of_processes;
//         int start = rpp * parId;
//         int end = start + rpp;
//         main_simpi->set_start(start);
//         main_simpi->set_end(end + tempForProcesses - 1) ;
//         if (Arow % number_of_processes != 0)
//         {

//             int leftover = Arow % number_of_processes;
//             if (parId < leftover)
//             {

//                 parId += (Arow - leftover);
//                 int start = parId;
//                 int end = start + 1;
//                 main_simpi->set_start(start);
//                 main_simpi->set_end(end + tempForProcesses - 1);
//                 for (int a = start; a < end; a++)
//                 {
//                     for (int b = 0; b < Arow; b++)
//                     {
//                         result->set(b + a * Arow,
//                         (other.get_algbera(b + a * Arow) + get_algbera(b + a * Arow)));
//                     }
//                 }
//             }
//         }
//         if (Acol != Brow)
//         {
//             printf("Error. Matrices can't be multiplied\n");
//         }
//         for (int a = start; a < end; a++)
//         {
//             for (int b = 0; b < Arow; b++)
//             {
//                 result->set(b + a * Arow,
//                         (other.get_algbera(b + a * Arow) + get_algbera(b + a * Arow)));
//             }
//         }
//         main_simpi->synch();
//         return *result;
//     }
// }
/*
Solves matrix subtraction in parallel and outputs the difference
matrix -> matrix
*/

matrix &matrix::subtract(matrix other)
{
    matrix *result = new matrix(xdim, other.get_y());
    int number_of_processes = main_simpi->get_synch_info()->par_count;
    int number_of_workstations = main_simpi->get_num_workstations();
    int tempForProcesses = number_of_processes;
    number_of_processes = number_of_processes * number_of_workstations;
    int parId = main_simpi->get_id();
    int workstationid = main_simpi->get_workstation_id() - 1 ;
    parId = workstationid * tempForProcesses + parId;
    printf("WORKSTATION ID = %d\n", workstationid);
    printf("parID ID = %d\n", parId);
    if (parId <= other.get_y())
    {
        if (number_of_processes > other.get_y())
        {
            number_of_processes = other.get_y();
        }
        int Arow = get_x();
        int Acol = get_y();
        int Brow = other.get_x();
        int Bcol = other.get_y();
        int rpp = Bcol / number_of_processes;
        int start = rpp * parId;
        int end = start + rpp;
        main_simpi->set_start(start);
        main_simpi->set_end(end + tempForProcesses - 1) ;
        if (Arow % number_of_processes != 0)
        {

            int leftover = Arow % number_of_processes;
            if (parId < leftover)
            {

                parId += (Arow - leftover);
                int start = parId;
                int end = start + 1;
                main_simpi->set_start(start);
                main_simpi->set_end(end + tempForProcesses - 1);
                for (int a = start; a < end; a++)
                {
                    for (int b = 0; b < Arow; b++)
                    {
                        result->set(b + a * Arow,
                        (other.get_algbera(b + a * Arow) - get_algbera(b + a * Arow)));
                    }
                }
            }
        }
        if (Acol != Brow)
        {
            printf("Error. Matrices can't be multiplied\n");
        }
        for (int a = start; a < end; a++)
        {
            for (int b = 0; b < Arow; b++)
            {
                result->set(b + a * Arow,
                        (other.get_algbera(b + a * Arow) - get_algbera(b + a * Arow)));
            }
        }
        main_simpi->synch();
        return *result;
    }
}

/*
Solves matrix scalar multiplication in parallel and outputs the resultant matrix
matrix -> matrix
*/
matrix &matrix::scalar_matrix_mult(int scaler)
{
    matrix *result = new matrix(get_x(), get_y());

    int number_of_processes = main_simpi->get_synch_info()->par_count;

    if (number_of_processes > get_x())
    {
        number_of_processes = get_x();
    }

    int rpp = get_y() / number_of_processes;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < get_y(); j++)
        {
            int pos = (get_y() * i + j);
            result->set(pos, get_algbera(pos) * scaler);
        }
    }
    return *result;
}

/*
Solves matrix equality and outputs the apt boolean value
matrix -> bool
*/
bool matrix::matrix_is_equal(matrix other)
{
    int number_of_processes = main_simpi->get_synch_info()->par_count;

    if (number_of_processes > get_x())
    {
        number_of_processes = get_x();
    }

    int Arow = get_x();
    int Acol = get_y();
    int Brow = other.get_x();
    int Bcol = other.get_y();
    int rpp = Arow / number_of_processes;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;

    if (Arow != Brow || Acol != Bcol)
    {
        return false;
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Acol; j++)
        {

            if (get_algbera(j + i * Acol) != other.get_algbera(j + i * Acol))
            {
                return false;
            }
        }
    }
    return true;
}

/*
Solves vector multiplication in parallel and outputs the resultant integer
vector -> int
*/

int vector::multiply(vector other)
{

    int sum = 0;
    if (get_size() != other.get_size())
    {
        printf("error: sizes don't match");
    }
    for (int i = 0; i < get_size(); i++)
    {
        sum += get(i) * other.get(i);
    }
    return sum;
}

/*
Solves vector scalar multiplication in parallel and outputs the resultant vector
vector,int -> vector
*/
vector &vector::scalar_vector_mult(int scaler)
{
    vector *result = new vector(get_size());
    int size = get_size();
    int rpp = size / main_simpi->get_synch_info()->par_count;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, get(i) * scaler);
    }
    return *result;
}

/*
Solves vector addition in parallel and outputs the resultant vector
vector-> vector
*/
vector &vector::add(vector other)
{
    vector *result = new vector(get_size());
    int size = get_size();
    int rpp = size / main_simpi->get_synch_info()->par_count;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, get(i) + other.get(i));
    }
    return *result;
}

/*
Solves vector subtraction in parallel and outputs the resultant vector
vector-> vector
*/
vector &vector::subtract(vector other)
{
    vector *result = new vector(get_size());
    int size = get_size();
    int rpp = size / main_simpi->get_synch_info()->par_count;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, get(i) - other.get(i));
    }
    return *result;
}

/*
Solves vector equality and outputs the apt boolean value
vector-> boolean
*/
bool vector::vector_is_equal(vector other)
{
    int size = get_size();
    int rpp = size / main_simpi->get_synch_info()->par_count;
    int start = rpp * main_simpi->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        if (get(i) != other.get(i))
        {
            return false;
        }
    }
    return true;
}

//Operator Overloading :-

//Multiplication
matrix &operator*(matrix &lhs, matrix &rhs)
{
    return lhs.multiply(rhs);
}
matrix &operator*(int lhs, matrix &rhs)
{
    return rhs.scalar_matrix_mult(lhs);
}
matrix &operator*(matrix &lhs, int rhs)
{
    return lhs.scalar_matrix_mult(rhs);
}
int operator*(vector &lhs, vector &rhs)
{
    return lhs.multiply(rhs);
}
vector &operator*(int lhs, vector &rhs)
{
    return rhs.scalar_vector_mult(lhs);
}
vector &operator*(vector &lhs, int rhs)
{
    return lhs.scalar_vector_mult(rhs);
}

//*=
void operator*=(matrix &lhs, matrix &rhs)
{
    lhs = lhs.multiply(rhs);
}
void operator*=(matrix &lhs, int rhs)
{
    lhs = lhs.scalar_matrix_mult(rhs);
}
void operator*=(vector &lhs, int rhs)
{
    lhs = lhs.scalar_vector_mult(rhs);
}

//Adition
matrix &operator+(matrix &lhs, matrix &rhs)
{
    return lhs.add(rhs);
}
vector &operator+(vector &lhs, vector &rhs)
{
    return lhs.add(rhs);
}

//Subtraction
matrix &operator-(matrix &lhs, matrix &rhs)
{
    return lhs.subtract(rhs);
}
vector &operator-(vector &lhs, vector &rhs)
{
    return lhs.subtract(rhs);
}

//+=
void operator+=(matrix &lhs, matrix &rhs)
{
    lhs = lhs.add(rhs);
}
void operator+=(vector &lhs, vector &rhs)
{
    lhs = lhs.add(rhs);
}

//-=
void operator-=(matrix &lhs, matrix &rhs)
{
    lhs = lhs.subtract(rhs);
}
void operator-=(vector &lhs, vector &rhs)
{
    lhs = lhs.subtract(rhs);
}

//Overloading ==
bool operator==(matrix &lhs, matrix &rhs)
{
    return lhs.matrix_is_equal(rhs);
}
bool operator==(vector &lhs, vector &rhs)
{
    return lhs.vector_is_equal(rhs);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

void pdDim_1_2(double *A, int dim, double *poly, int *lambdas)
{
    /*
 dim == 1 only when polyDet is called from the eigen value method,
 this will most likely not exist once the eigen value method is made and
 1x1 matrices are easily calculated on their own.
   */
    if (dim == 1)
    {
        poly[0] = A[0];
        poly[1] = -1;
        return;
    }
    // base case
    else
    {
        double polyad[3];
        double polybc[3];
        //mulitply ad
        if (lambdas[0] == 0 && lambdas[3] == 0)
        {
            polyad[0] = A[0] * A[3];
            polyad[1] = 0;
            polyad[2] = 0;
        }
        else if (lambdas[0] == 1 && lambdas[3] == 0)
        {
            polyad[0] = A[0] * A[3];
            if (A[3] == 0)
            {
                polyad[1] = 0;
            }
            else
            {
                polyad[1] = -1 * A[3];
            }
            polyad[2] = 0;
        }
        else if (lambdas[0] == 0 && lambdas[3] == 1)
        {
            polyad[0] = A[0] * A[3];
            if (A[0] == 0)
            {
                polyad[1] = 0;
            }
            else
            {
                polyad[1] = -1 * A[0];
            }
            polyad[2] = 0;
        }
        else
        {
            polyad[0] = A[0] * A[3];
            if (A[0] == 0 && A[3] == 0)
            {
                polyad[1] = 0;
            }
            else
            {
                polyad[1] = -1 * (A[0] + A[3]);
            }
            polyad[2] = 1;
        }
        //multiply bc
        if (lambdas[1] == 0 && lambdas[2] == 0)
        {
            polybc[0] = A[1] * A[2];
            polybc[1] = 0;
            polybc[2] = 0;
        }
        else if (lambdas[1] == 1 && lambdas[2] == 0)
        {
            polybc[0] = A[1] * A[2];
            if (A[2] == 0)
            {
                polybc[1] = 0;
            }
            else
            {
                polybc[1] = -1 * A[2];
            }
            polybc[2] = 0;
        }
        else if (lambdas[1] == 0 && lambdas[2] == 1)
        {
            polybc[0] = A[1] * A[2];
            if (A[1] == 0)
            {
                polybc[1] = 0;
            }
            else
            {
                polybc[1] = -1 * A[1];
            }
            polybc[2] = 0;
        }
        else
        {
            polybc[0] = A[1] * A[2];
            if (A[1] == 0 && A[2] == 0)
            {
                poly[1] = 0;
            }
            else
            {
                polybc[1] = -1 * (A[1] + A[2]);
            }
            polybc[2] = 1;
        }
        // subtracting ad-bc and putting result in poly
        if (polybc[0] == 0)
        {
            poly[0] = polyad[0];
        }
        else
        {
            poly[0] = polyad[0] - polybc[0];
        }
        if (polybc[1] == 0)
        {
            poly[1] = polyad[1];
        }
        else
        {
            poly[1] = polyad[1] - polybc[1];
        }
        if (polybc[2] == 0)
        {
            poly[2] = polyad[2];
        }
        else
        {
            poly[2] = polyad[2] - polybc[2];
        }
        return;
    }
}

void polyDet_aux(double *A, int dim, double *poly, int *lambdas)
{
    if (dim == 1 || dim == 2)
    {
        pdDim_1_2(A, dim, poly, lambdas);
        return;
    }
    int sign = 0; //0 to add 1 to subtract
    // for each value in the top row
    for (int x = 0; x < dim; x++)
    {
        // make new matrix
        double subA[(dim - 1) * (dim - 1)];
        int idx = 0;
        int newLambdas[(dim - 1) * (dim - 1)];
        for (int i = 1; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                if (j != x)
                {
                    subA[idx] = A[j + dim * i];
                    if (lambdas[j + dim * i] == 1)
                    {
                        newLambdas[idx] = 1;
                    }
                    else
                    {
                        newLambdas[idx] = 0;
                    }
                    idx++;
                }
            }
        }
        double subPoly[dim + 1];
        for (int sp = 0; sp < dim + 1; sp++)
        {
            subPoly[sp] = 0;
        }
        polyDet_aux(subA, dim - 1, subPoly, newLambdas);
        // add and subtract sub-determinants
        if (sign == 0)
        {
            for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
            {
                poly[i] += subPoly[i] * A[x];
            }
            if (lambdas[x] == 1)
            {
                for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
                {
                    poly[i + 1] -= subPoly[i];
                }
            }
            else
            {
                // highest possible power is of lambda is len(poly)-2
                poly[dim + 1] = 0;
            }
            sign = 1;
        }
        else
        {
            double subtractPoly[dim + 1];
            for (int i = 0; i < dim + 1; i++)
            {
                subtractPoly[i] = 0;
            }
            for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
            {
                subtractPoly[i] += subPoly[i] * A[x];
            }
            if (lambdas[x] == 1)
            {
                for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
                {
                    subtractPoly[i + 1] -= subPoly[i];
                }
            }
            else
            {
                // highest possible power is of lambda is len(poly)-2
                subtractPoly[dim + 1] = 0;
            }
            for (int i = 0; i < (sizeof(subPoly) / sizeof(double)) + 1; i++)
            {
                poly[i] -= subtractPoly[i];
            }
            sign = 0;
        }
    }
}

void pre_polyDet_aux(double *A, int dim, double *poly, int *lambdas, int sign, int x)
{
    double subA[(dim - 1) * (dim - 1)];
    int idx = 0;
    int newLambdas[(dim - 1) * (dim - 1)];
    for (int i = 1; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (j != x)
            {
                subA[idx] = A[j + dim * i];
                if (lambdas[i + dim * j] == 1)
                {
                    newLambdas[idx] = 1;
                }
                else
                {
                    newLambdas[idx] = 0;
                }
                idx++;
            }
        }
    }
    double subPoly[dim + 1];
    for (int sp = 0; sp < dim + 1; sp++)
    {
        subPoly[sp] = 0;
    }
    polyDet_aux(subA, dim - 1, subPoly, newLambdas);
    // add and subtract sub-determinants
    if (sign == 0)
    {
        for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
        {
            poly[i] += subPoly[i] * A[x];
        }
        if (lambdas[x] == 1)
        {
            for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
            {
                poly[i + 1] -= subPoly[i];
            }
        }
        else
        {
            // highest possible power is of lambda is len(poly)-2
            poly[dim + 1] = 0;
        }
    }
    else
    {
        double subtractPoly[dim + 1];
        for (int i = 0; i < dim + 1; i++)
        {
            subtractPoly[i] = 0;
        }
        for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
        {
            subtractPoly[i] += subPoly[i] * A[x];
        }
        if (lambdas[x] == 1)
        {
            for (int i = 0; i < sizeof(subPoly) / sizeof(double); i++)
            {
                subtractPoly[i + 1] -= subPoly[i];
            }
        }
        else
        {
            // highest possible power is of lambda is len(poly)-2
            subtractPoly[dim + 1] = 0;
        }
        for (int i = 0; i < (sizeof(subPoly) / sizeof(double)) + 1; i++)
        {
            poly[i] -= subtractPoly[i];
        }
    }
}

double power(double num, int pow)
{
    double temp = num;
    for (int i = 0; i < pow - 1; i++)
    {
        temp = temp * num;
    }
    return temp;
}
/*
The input poly is a float array with each index representing a multiple of that
index-power of x in the polynomial equation d+cx+bx^2+..+ax^(len(P)-1)

The method polyRoots solve for the roots of this equation and returns a float
array of all zeroes of P

This method is based off of The Rational Zero Theorem
*/
void polyRoots(double *poly, double *roots, int polySize)
{
    /*
      number roots cant exceed number of factors of the constant term of the polynomial which is always less than
      or equal to the number itself
      the rest will be filled in with 0s as 0 cannot be a valid eigan value

      P = factors of the constant term
      Q = factors of the leading coefficient
    */
    // round decimal place
    for (int i = 0; i < polySize; i++)
    {
        poly[i] = floor(poly[i] * 10000) / 10000;
    }
    if (poly[0] == 0)
    {
        // take out a factor of lambda and find the new solution with a root of zero
        double subpoly[polySize - 1];
        for (int i = 1; i < polySize; i++)
        {
            subpoly[i - 1] = poly[i];
        }
        double subroots[polySize - 1];
        polyRoots(subpoly, subroots, polySize - 1);
        roots[0] = 0;
        for (int i = 1; i < polySize; i++)
        {
            roots[i] = subroots[i - 1];
        }
        return;
    }
    int P = (int)abs(poly[0]);
    int Q = (int)abs(poly[polySize - 1]);
    // decimal adjustment
    int factor = 0;
    while (abs(P - abs(poly[0])) > 0.00001 || abs(Q - abs(poly[polySize - 1])) > 0.00001)
    {
        for (int i = 0; i < polySize; i++)
        {
            poly[i] = 10 * poly[i];
        }
        factor++;
        P = (int)abs(poly[0]);
        Q = (int)abs(poly[polySize - 1]);
    }
    /*
     these max's are causing a limit on how many roots can be found
     but they are causing an error when creating the factorsofX arrays
     when they are too big
    */
    int maxnumPs = P / 2;
    int maxnumQs = Q / 2;
    // find factors of P and Q
    int factorsOfP[maxnumPs];
    int fopIDX = 0;
    for (int i = 1; i < P + 1; i++)
    {
        if (P % i == 0)
        {
            factorsOfP[fopIDX] = i;
            fopIDX++;
        }
    }
    int factorsOfQ[maxnumQs];
    int foqIDX = 0;
    for (int i = 1; i < Q + 1; i++)
    {
        if (Q % i == 0)
        {
            factorsOfQ[foqIDX] = i;
            foqIDX++;
        }
    }
    // create our P/Q values
    int pos_rootsIDX = 0;
    double possible_roots[foqIDX * fopIDX * 2];
    for (int i = 0; i < foqIDX; i++)
    {
        for (int j = 0; j < fopIDX; j++)
        {
            possible_roots[pos_rootsIDX] = (double)factorsOfP[j] / factorsOfQ[i];
            pos_rootsIDX++;
            possible_roots[pos_rootsIDX] = -1 * (double)factorsOfP[j] / factorsOfQ[i];
            pos_rootsIDX++;
        }
    }
    // restore poly if necessary---testing
    if (factor != 0)
    {
        for (int i = 0; i < factor; i++)
        {
            for (int j = 0; j < polySize; j++)
            {
                poly[j] = poly[j] / 10;
            }
        }
    }
    // evaluate possible roots
    int rootsIDX = 0;
    for (int i = 0; i < pos_rootsIDX; i++)
    {
        double evaluation;
        for (int j = 0; j < polySize; j++)
        {
            if (j == 0)
            {
                evaluation = poly[0];
            }
            else
            {
                evaluation += (double)poly[j] * power(possible_roots[i], j);
            }
        }
        int already_root = 0;
        for (int j = 0; j < rootsIDX; j++)
        {
            if (roots[j] == possible_roots[i])
            {
                already_root = 1;
            }
        }
        if (abs(evaluation) < 0.00001 && already_root == 0)
        {

            roots[rootsIDX] = possible_roots[i];
            rootsIDX++;
        }
    }
}

/*
spit up work and solve for for the polynomial given from the determinent of
the equation A-lambda*I
*/
void polyDet(int polyFD, int par_id, int par_count, int dim, double *A, double *poly, int *lambdas)
{
    if (dim < 3)
    {
        pdDim_1_2(A, dim, poly, lambdas);
        return;
    }
    int work = dim / par_count;
    if (work < 1)
    {
        work = 1;
    }
    int start = par_id * work;
    int end = start + work;
    int xtra = dim - (work * par_count);
    if (xtra != 0 && par_count - xtra - 1 < par_id)
    {
        start += (par_id - (par_count - xtra));
        end = start + work + 1;
    }
    // split up the remaining processes
    for (int x = start; x < end; x++)
    {
        int sign = 0;
        if (x % 2 != 0)
        {
            sign = 1;
        }
        pre_polyDet_aux(A, dim, poly, lambdas, sign, x);
    }
}
double *eigenvalue(matrix Anaught, int par_id, int par_count)
{
    int dim = Anaught.get_x();
    if (Anaught.get_x() != Anaught.get_y())
    {
        throw std::invalid_argument("matrix must be a square matrix\n");
    }
    // create lambda matrix and copy of A
    int lambdas[dim * dim];
    matrix A = matrix(dim, dim);
    for (int i = 0; i < dim * dim; i++)
    {
        A.arr[i] = Anaught.arr[i];
        if (i % (dim + 1))
        {
            lambdas[i] = 0;
        }
        else
        {
            lambdas[i] = 1;
        }
    }
    // init poly memory
    double *poly;
    double *roots;
    int polyFD;
    int rootsFD;
    if (par_id == 0)
    {
        polyFD = shm_open("eig_val_poly", O_RDWR | O_CREAT, 0777);
        if (polyFD == -1)
        {
            printf("shared memory failed in eigen value calculation\n");
        }
        ftruncate(polyFD, sizeof(double) * (dim + 1));
        poly = (double *)mmap(NULL, sizeof(double) * (dim + 1), PROT_READ | PROT_WRITE, MAP_SHARED, polyFD, 0);
        rootsFD = shm_open("eig_val_roots", O_RDWR | O_CREAT, 0777);
        if (rootsFD == -1)
        {
            printf("shared memory failed in eigen value calculation\n");
        }
        ftruncate(rootsFD, sizeof(double) * (dim + 1));
        roots = (double *)mmap(NULL, sizeof(double) * (dim + 1), PROT_READ | PROT_WRITE, MAP_SHARED, rootsFD, 0);
        for (int i = 0; i < dim + 1; i++)
        {
            poly[i] = 0;
            roots[i] = 0;
        }
    }
    main_simpi->synch();
    if (par_id != 0)
    {
        // sleep(2);
        polyFD = shm_open("eig_val_poly", O_RDWR, 0777);
        poly = (double *)mmap(NULL, sizeof(double) * (dim + 1), PROT_READ | PROT_WRITE, MAP_SHARED, polyFD, 0);
        rootsFD = shm_open("eig_val_roots", O_RDWR, 0777);
        roots = (double *)mmap(NULL, sizeof(double) * (dim + 1), PROT_READ | PROT_WRITE, MAP_SHARED, rootsFD, 0);
    }
    // get the polynomial from the determinent
    if (dim > par_id)
    {
        polyDet(polyFD, par_id, par_count, dim, A.arr, poly, lambdas);
    }
    main_simpi->synch();
    if (par_id == 0)
    {
        polyRoots(poly, roots, dim + 1);
    }
    shm_unlink("eig_val_poly");
    main_simpi->synch();
    return roots;
}
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <string.h>

/*
 takes a matrix in row echelon form and solves the system of equations set
 equal to zero
 */
void solve_REF_equations(double *M, double *answer, int dim)
{
    int lastrc = dim - 1;
    if (M[lastrc + dim * lastrc] != 0)
    {
        answer[lastrc] = 0;
    }
    else
    {
        answer[lastrc] = 1;
    }
    for (int row = lastrc - 1; row >= 0; row--)
    {
        int new_col = row;
        double known_sum = 0;
        for (int i = new_col + 1; i < dim; i++)
        {
            known_sum += (M[i + row * dim] * answer[i]);
        }
        if (known_sum == 0 || M[new_col + row * dim] == 0)
        {
            answer[new_col] = 0;
        }
        else
        {
            answer[new_col] = (-1 * known_sum) / M[new_col + row * dim];
        }
    }
}

// helper method of rowEclon()
void swapRows(double *M, int i, int j, int dim)
{
    for (int k = 0; k < dim; k++)
    {
        double temp = M[k + i * dim];
        M[k + i * dim] = M[k + j * dim];
        M[k + j * dim] = temp;
    }
}

// takes a double array matrix and puts it into row echelon form
void rowEchelon(double *M, int dim)
{
    for (int k = 0; k < dim; k++)
    {
        int i_max = k;
        int v_max = M[k + i_max * dim];
        for (int i = k + 1; i < dim; i++)
        {
            if (abs(M[k + i * dim]) > v_max)
            {
                v_max = M[k + i * dim];
                i_max = i;
            }
        }
        if (i_max != k)
        {
            swapRows(M, k, i_max, dim);
        }
        for (int i = k + 1; i < dim; i++)
        {
            double f = M[k + i * dim] / M[k + k * dim];
            for (int j = k + 1; j < dim; j++)
            {
                M[j + i * dim] -= M[j + k * dim] * f;
                M[k + i * dim] = 0;
            }
        }
    }
    for (int i = 0; i < dim * dim; i++)
    {
        if (abs(M[i]) < 0.00001)
        {
            M[i] = 0;
        }
    }
}

/*
   matrix A: the matrix that you want the eigen values for
   double *V: the double array you want the vectors to end up in
   par_id: the process id if each individual processes
   par_count: the total number of processes running
*/
void eigenvector(matrix A, double *V, int par_id, int par_count)
{
    // check valid matrix
    if (A.get_x() != A.get_y())
    {
        throw std::invalid_argument("matrix must be a square matrix");
    }
    int dim = A.get_x();
    matrix *copyA = new matrix(dim, dim);
    for (int i = 0; i < dim * dim; i++)
    {
        copyA->arr[i] = A.arr[i];
    }
    double *eigval = eigenvalue(*copyA, par_id, par_count);
    if (par_id == 0)
    {
        for (int i = 0; i < dim; i++)
        {
            double AminL[dim * dim];
            for (int j = 0; j < dim; j++)
            {
                for (int k = 0; k < dim; k++)
                {
                    if (j == k)
                    {
                        AminL[k + dim * j] = A.arr[k + dim * j] - eigval[i];
                    }
                    else
                    {
                        AminL[k + j * dim] = A.arr[k + j * dim];
                    }
                }
            }
            rowEchelon(AminL, dim);
            double evect[dim];
            solve_REF_equations(AminL, evect, dim);
            for (int j = 0; j < dim; j++)
            {
                V[j + i * dim] = evect[j];
            }
        }
    }
}
