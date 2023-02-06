#include "def.h"



int CRACKTHRESHOLD = 0, biggestWordSize = 0;
Node *root;
float **points, *distances, **queries;
int *ids, *gc;
vector<string> query_words;
char **words;

#ifdef L2
    #define distance(point,query,dimension) (tools::L2Distance(point,query,dimension))
#elif L1
     #define distance(point,query,dimension) (tools::L1Distance(point,query,dimension))
#elif ED
     #define distance(point,query,dimension) (tools::EditDistance(point,query))
#endif


namespace tools {

    void loadData(const char *filename, float ***centers, int *numObjects, int *dimension)
    {
        int dim, num_obj, func, id=0;
        float x;

        ifstream inp(filename);
        inp >> dim >> num_obj >> func;
        *numObjects = num_obj;
        *dimension = dim;
        cout << dim << " " << num_obj << " " << func << endl;

        *centers = (float **)malloc(*numObjects*sizeof(float *));
        for(int i=0;i<*numObjects;i++)
        {
            (*centers)[i] = (float *)malloc(dim*sizeof(float));
        }

        int i = 0;
        string line;
        getline(inp, line);
        while(getline(inp, line)) { 
            stringstream ss(line);
            int j = 0;
            while(ss && !ss.eof()) {
                ss >> x;
                (*centers)[i][j] = x;
                j++;
            }
            i++;
        }
        inp.close();
    }

    void loadQueries(const char *filename, float ***queries, float **radius,int *numQueries,int dim)
    {
        int num_qry=0, func, id=0;
        float x;

        ifstream inp(filename);
        inp >> num_qry;  
        *numQueries = num_qry;

        *radius = (float *)malloc(*numQueries*sizeof(float));	
        *queries = (float **)malloc(*numQueries*sizeof(float *));
        for(int i=0;i<*numQueries;i++)
        {
            (*queries)[i] = (float *)malloc(dim*sizeof(float));
        }

        int i = 0;
        string line;
        getline(inp, line);
        while(getline(inp, line)) { 
            stringstream ss(line);
            int j = 0;
            while(ss && !ss.eof()) {
                ss >> x;
                if(j==0)
                {
                    (*radius)[i] = x;
                }
                else if(j<=dim)
                {
                    (*queries)[i][j-1] = x;
                }
                j++;
            }
            i++;
        }
        inp.close();
    }

    char** loadDataString(const char *filename, int *numObjects)
    {
        char *ptr,*top;
        FILE *f;
        struct stat sdata;
        unsigned long dn;
        int DBlen;
        char *pals;    /* words all together */
        char **ptrs;   /* pointers to each word */

        f = fopen (filename,"r");
        stat (filename,&sdata);
        pals = (char *)malloc (sdata.st_size*sizeof(char));
        DBlen = sdata.st_size;
        fread (pals,sdata.st_size,1,f);
        fclose (f);
        ptr = pals; top = ptr + sdata.st_size;
        dn = 0;
        while (ptr < top) 
        { 
            while (*ptr != '\n') ptr++;
            dn++; *ptr++ = 0;
        }
        ptrs = (char**)malloc ((dn+1)*sizeof(char*));
        dn = 0; ptr = pals;
        ptrs[0] = NULL;
        while (ptr < top) 
        { 
            if(strlen(ptr) > biggestWordSize)
                biggestWordSize = strlen(ptr);
            ptrs[dn] = ptr;
            dn++;
            while (*ptr++);
        }
        *numObjects = dn;
        gc = (int*)realloc(gc,(biggestWordSize+1)*sizeof(int)); 
        return ptrs;
    }

    void loadQueriesString(const char *filename,int *numQueries, float **radius)
    {
        int num_qry=0, func, id=0;
        char* x , *word;
        int counter = 0;

        ifstream inp(filename);
        inp >> num_qry;  
        *numQueries = num_qry;

        *radius = (float *)malloc(*numQueries*sizeof(float));	
        query_words.resize(num_qry);
        
        string line;
        int i = 0;
        while(std::getline(inp, line))
        {
            string str1, str2;
            stringstream s(line);

            s>>str1>>str2;
            if(str1.size() >0)
            {
                (*radius)[i] = stof(str1);
                query_words[i]=str2;
                i++;
            }
        }
        inp.close();
    }

    static float L2Distance (float *p1, float *p2, int k)
    { 
        register int i;
        float tot = 0,dif;

        for (i=0;i<k;i++) 
        { 
            dif = (p1[i]-p2[i]);
            tot += dif*dif;
        }
        return sqrt(tot);
    }

    static float L1Distance (float *p1, float *p2, int k)
    { 
        register int i;
        float tot = 0,dif;
        for (i=0;i<k;i++) 
        { 
            dif = (p1[i]-p2[i]);
            if (dif < 0) dif = -dif;
            tot += dif;
        }
        return tot;
    }

    static int EditDistance (char *p1, char *p2)
    { 
        register int pc,nc,j,i; 
        register char cc; 
        register int m = strlen(p1);
        int *c;
        if (biggestWordSize < m) 
	    {  
            biggestWordSize = m; 
	        gc = (int*)realloc(gc,(m+1)*sizeof(int)); 
	    }
        c = gc;
        nc = m;
        p1--;
        for (j=0;j<=m;j++) c[j] = j;
        for (i=0;(cc=p2[i]);i++)
        { 
            pc = i; nc = i+1; 
            for (j=1;j<=m;j++) 
            { 
                if (c[j] < nc) nc = c[j];
                pc += (cc != p1[j]);
                if (pc <= nc) nc = pc; else nc++; 
                pc = c[j]; c[j] = nc; 
            } 
        } 
        return nc;
    }

}

#ifdef ED
    namespace standard{

        int crackInTwo(char **points, float *distances, int start, int end, float epsilon, char* query, int *found_count, int *ids, float *maxDistance)
        {
            int i = start, j = end;
            while(true)
            {
                
                while(distances[i] <= epsilon && i < end+1){
                    if(distances[i] > *maxDistance)
                        (*maxDistance) = distances[i];
                    i++;
                    (*found_count)++;
                }
                while(distances[j] > epsilon && j>start-1){
                    if(distances[j] > *maxDistance)
                        (*maxDistance) = distances[j];
                    j--;   
                }
                if(i>=j)
                    break;
                swap(points[i], points[j]);       
                swap(distances[i], distances[j]);
                swap(ids[i], ids[j]);
            }
            return j;
        }

        void searchAndCrack(Node *T, char *query, float epsilon, char **points, float **distances, int *result, int *ids, int dimension, int queryId)
        {   
            int i, j, crack;
            float maxDistance, dist;

            if(T->right == NULL  && T->left == NULL)
            {
                    
                for(i = T->start; i <= T->end; i++)
                    (*distances)[i] = distance(points[i],query,0);

                maxDistance = -1.0;
                crack = crackInTwo(points,*distances,T->start,T->end,epsilon, query, result, ids, &maxDistance);

                T->rec = query;
                T->epsilon = epsilon;
                if(crack >= T->start && T->end >= crack+1)
                {
                    T = T->insertLeftString(T, NULL, T->start, crack,queryId);
                    T = T->insertRightString(T, NULL, crack+1, T->end,queryId );
                }
            }
            else
            {
                dist = distance(T->rec,query,0);
                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++)
                        {
                            (*result)++;
                        }
                        searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                    }
                    else
                    {
                        searchAndCrack(T->left, query, epsilon, points, distances, result, ids, dimension, queryId);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                        }
                    }
                }
                else
                {
                    searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                }
            }
        }
    }

    namespace standardMediocre
    {

        int crackInTwoMediocre(char **points, float *distances, int start, int end, float epsilon, float newE, char* query, int *found_count, int *ids, float *maxDistance)
        {
            int i = start, j = end;
            while(true)
            {
                while(distances[i] <= newE && i < end+1){
                    if(distances[i] <= epsilon)
                    {
                        if(distances[i] > *maxDistance)
                            (*maxDistance) = distances[i];
                        (*found_count)++;
                    }
                    i++;
                }

                while(distances[j] > newE && j>start-1){
                    if(distances[j] <= epsilon)
                    {
                        if(distances[j] > *maxDistance)
                            (*maxDistance) = distances[j];
                        (*found_count)++;
                    } 
                    j--;   
                }
                if(i>=j)
                    break;
                swap(points[i], points[j]);       
                swap(distances[i], distances[j]);
                swap(ids[i], ids[j]);
            }

            return j;
        }

        void searchAndCrackMediocre(Node *T, char *query, float epsilon, char **words, float **distances, int *result, int *ids, int dimension, int queryId)
        {   
            int i, j, crack;
            float maxDistance, dist, median;
            vector<int> rnd_dists;

            if(T->right == NULL  && T->left == NULL) 
            {
                for(i = T->start; i <= T->end; i++)
                    (*distances)[i] = distance(words[i],query,0);

                rnd_dists.clear();
                rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));

                median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));
                maxDistance = -1.0;
                crack = crackInTwoMediocre(words, *distances, T->start,T->end,epsilon, median, query, result, ids, &maxDistance);

                T->rec = query;
                T->epsilon = median;
                
                if(crack >= T->start && T->end >= crack+1)
                {
                    T = T->insertLeftString(T, NULL, T->start, crack,queryId);
                    T = T->insertRightString(T, NULL, crack+1, T->end,queryId );
                }
            }
            else
            {
                dist = distance(T->rec,query,0);
                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(int j = T->left->start; j <= T->left->end; j++) // TODO remove for and just add
                        {
                            (*result)++;
                        }
                        searchAndCrackMediocre(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                    }
                    else
                    {
                        searchAndCrackMediocre(T->left, query, epsilon, words, distances, result, ids, dimension, queryId);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrackMediocre(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                        }
                    }
                }
                else
                {
                    searchAndCrackMediocre(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                }
            }
        }

        void searchAndCrackMediocreThreshold(Node *T, char *query, float epsilon, char **words, float **distances, int *result, int *ids, int dimension, int queryId)
        {   
            int i, j, crack;
            float maxDistance, dist, median;
            vector<int> rnd_dists;

            if(T->right == NULL  && T->left == NULL)
            {
                for(i = T->start; i <= T->end; i++)
                    (*distances)[i] = distance(words[i],query,0);

                if((T->end - T->start + 1) <= CRACKTHRESHOLD){
                    for(int i = T->start; i <= T->end ; i++){
                        if((*distances)[i] <= epsilon){
                            (*result)++;
                        }
                    }
                }
                else
                {
                    rnd_dists.clear();
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));

                    median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));

                    maxDistance = -1.0;
                    crack = crackInTwoMediocre(words, *distances, T->start,T->end,epsilon, median, query, result, ids, &maxDistance);

                    T->rec = query;
                    T->epsilon = median;
                    
                    if(crack >= T->start && T->end >= crack+1)
                    {
                        T = T->insertLeftString(T, NULL, T->start, crack,queryId);
                        T = T->insertRightString(T, NULL, crack+1, T->end,queryId );
                    }
                }
            }
            else
            {
                dist = distance(T->rec,query,0);
                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++)
                        {
                            (*result)++;
                        }
                        searchAndCrackMediocreThreshold(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                    }
                    else
                    {
                        searchAndCrackMediocreThreshold(T->left, query, epsilon, words, distances, result, ids, dimension, queryId);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrackMediocreThreshold(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                        }
                    }
                }
                else
                {
                    searchAndCrackMediocreThreshold(T->right, query, epsilon, words, distances, result, ids, dimension, queryId);
                }
            }
        }

    }

    namespace standard_cache_sort
    {

        int partition(int low, int high, int* lp) 
        { 
            if (distances[low] > distances[high]) 
            {
                swap(ids[low], ids[high]); 
                swap(distances[low], distances[high]); 
                swap(words[low], words[high]); 
            }
            int j = low + 1; 
            int g = high - 1, k = low + 1;
            float p = distances[low], q = distances[high]; 
            while (k <= g) { 
        
                if (distances[k] < p) { 
                    swap(ids[k], ids[j]); 
                    swap(distances[k], distances[j]); 
                    swap(words[k], words[j]); 
                    j++; 
                } 
                else if (distances[k] >= q) { 
                    while (distances[g] > q && k < g) 
                        g--; 
                    swap(ids[k], ids[g]); 
                    swap(distances[k], distances[g]); 
                    swap(words[k], words[g]); 
                    g--; 
                    if (distances[k] < p) { 
                        swap(ids[k], ids[j]);
                        swap(distances[k], distances[j]);
                        swap(words[k], words[j]);
                        j++; 
                    } 
                } 
                k++; 
            } 
            j--; 
            g++; 
        
            // bring pivots to their appropriate positions. 
            swap(ids[low], ids[j]); 
            swap(ids[high], ids[g]); 

            swap(distances[low], distances[j]); 
            swap(distances[high], distances[g]);

            swap(words[low], words[j]); 
            swap(words[high], words[g]);

            // returning the indices of the pivots. 
            *lp = j; // because we cannot return two elements 
            // from a function. 
        
            return g; 
        } 

        void DualPivotQuickSort(int low, int high) 
        { 
            if (low < high) { 
                int lp, rp; 
                rp = partition(low, high, &lp);
                DualPivotQuickSort(low, lp - 1); 
                DualPivotQuickSort(lp + 1, rp - 1); 
                DualPivotQuickSort(rp + 1, high); 
            } 
        } 

        int binarySearchFirstBigger(float **R, int start, int end, float key)
        {
            if(start >= end)
                return 0;
            
            int low = start;
            int high = end;
            
            int ans = end;
            bool flag = false;
            while (low <= high) {
                
                int mid = (high+low)/2;
                double midVal = (*R)[mid];
                
                if (midVal < key) {
                    low = mid + 1;
                }
                else if (midVal >= key) {
                    ans = mid;
                    high = mid - 1;
                }
            }
            if(low > end)
                return end+1;
            
            return ans;
        }

        int binarySearchLastSmaller(float **R, int start, int end, float key)
        {
            if(start >= end)
                return 0;
            
            int low = start;
            int high = end;
            
            int ans = end;
            bool flag = false;
            while (low <= high) {
                
                int mid = (high+low)/2;
                double midVal = (*R)[mid];
                
                if (midVal <= key) {
                    low = mid + 1;
                    ans = mid;
                }
                else if (midVal > key) {
                    high = mid - 1;
                }
            }
            if(high < start)
                return start-1;
            return ans;
        }
        
        
        void searchAndCrackMediocreThresholdSortBinaryOpt(Node *T, char *query, float epsilon, char **words, float **distances, int *result, int *ids, int dimension, int queryId, Node* prev, float qDist)
        {   
            int i, j, crack;
            float maxDistance, dist, median;
            vector<int> rnd_dists;
            int low, high;

            if(T->right == NULL  && T->left == NULL)
            {
                if((T->end - T->start + 1) <= CRACKTHRESHOLD)
                {
                    float qDist = distance(prev->rec, query,0);
                    if(T->leftchild == 1)
                    {
                        if(qDist > epsilon)                                                                 // CASE L1
                        {
                            if((*distances)[T->start] >= qDist - epsilon)
                                low = T->start;
                            else
                                low = standard_cache_sort::binarySearchFirstBigger(distances, T->start, T->end, qDist-epsilon);
                            if(qDist >= prev->epsilon)  // POTENTIAL BUG!!!!!     
                                high = T->end;
                            else
                            {
                                high = standard_cache_sort::binarySearchLastSmaller(distances, T->start, T->end, qDist + epsilon);
                            }

                            for(i = low; i <= high; i++)
                            {
                                if( distance(words[i],query,0) <= epsilon)
                                {
                                    (*result)++;
                                }
                            }
                        }
                        else            // CASE L2
                        {
                            if(qDist + prev->epsilon <= epsilon)
                                low = T->end+1;
                            else
                                low = standard_cache_sort::binarySearchFirstBigger(distances, T->start, T->end, epsilon-qDist);

                            for(i = T->start; i < low; i++)
                            {
                                (*result)++;
                            }
                            high = T->end;

                            for(i = low; i <= T->end; i++)
                            {
                                if(distance(words[i],query,0) <= epsilon)
                                {
                                    (*result)++;
                                }
                            }
                        }
                    } 
                    else
                    {
                        if(qDist > epsilon+prev->epsilon)
                        {
                            // cout << "CASE R1\n";
                            if((*distances)[T->start] >= qDist + epsilon)
                                low = T->start;
                            else
                                low = standard_cache_sort::binarySearchFirstBigger(distances, T->start, T->end, qDist-epsilon);
                        }
                        else
                        {
                            // cout << "CASE R2\n";
                            low = T->start;
                        }

                        if((*distances)[T->end] <= qDist + epsilon)
                            high = T->end;
                        else
                            high = standard_cache_sort::binarySearchLastSmaller(distances, T->start, T->end, qDist+epsilon);

                        for(i = low; i <= high; i++)
                        {
                            if(distance(words[i],query,0) <= epsilon)
                            {
                                (*result)++;
                            }  
                        }

                    }
                }
                else{
                    for(i = T->start; i <= T->end; i++)
                        (*distances)[i] = distance(words[i],query,0);

                    rnd_dists.clear();
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));

                    median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));
                    
                    maxDistance = -1.0;
                    crack = standardMediocre::crackInTwoMediocre(words, *distances, T->start,T->end,epsilon, median, query, result, ids, &maxDistance);

                    T->rec = query;
                    T->epsilon = median;


                    if(crack >= T->start && T->end >= crack+1)
                    {
                        T->id = queryId;
                        T = T->insertLeftString(T, NULL, T->start, crack,queryId, 1);
                        T = T->insertRightString(T, NULL, crack+1, T->end,queryId, 0 );
                        if(crack - T->start + 1 <= CRACKTHRESHOLD)
                        {
                            standard_cache_sort::DualPivotQuickSort(T->start, crack );          // Sort piece in ascending order
                        }

                        if(T->end - crack <= CRACKTHRESHOLD)
                        {
                            standard_cache_sort::DualPivotQuickSort(crack+1, T->end);  // Sort piece in ascending order
                        }

                    }
                }
            }
            else
            {
                dist = distance(T->rec,query,0);

                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++) 
                        {
                            (*result)++;
                        }
                        searchAndCrackMediocreThresholdSortBinaryOpt(T->right, query, epsilon, words, distances, result, ids, dimension, queryId, T, dist);
                    }
                    else
                    {
                        searchAndCrackMediocreThresholdSortBinaryOpt(T->left, query, epsilon, words, distances, result, ids, dimension, queryId, T, dist);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrackMediocreThresholdSortBinaryOpt(T->right, query, epsilon, words, distances, result, ids, dimension, queryId, T, dist);
                        }
                    }
                }
                else
                {
                    searchAndCrackMediocreThresholdSortBinaryOpt(T->right, query, epsilon, words, distances, result, ids, dimension, queryId, T, dist);
                }
            }
        }

    }

#else

    namespace standard{

        int crackInTwo(float **points, float *distances, int start, int end, float epsilon, float* query, int *found_count, int *ids, float *maxDistance)
        {
            int i = start, j = end;
            while(true)
            {
                
                while(distances[i] <= epsilon && i < end+1){
                    if(distances[i] > *maxDistance)
                        (*maxDistance) = distances[i];
                    i++;
                    (*found_count)++;
                }
                while(distances[j] > epsilon && j>start-1){
                    if(distances[j] > *maxDistance)
                        (*maxDistance) = distances[j];
                    j--;   
                }
                if(i>=j)
                    break;
                swap(points[i], points[j]);       
                swap(distances[i], distances[j]);
                swap(ids[i], ids[j]);
            }
            return j;
        }

        void searchAndCrack(Node *T, float *query, float epsilon, float **points, float **distances, int *result, int *ids, int dimension, int queryId)
        {   
            int i, j, crack;
            float maxDistance, dist;

            if(T->right == NULL  && T->left == NULL)
            {
                    
                for(i = T->start; i <= T->end; i++)
                    (*distances)[i] = distance(points[i],query,dimension);

                maxDistance = -1.0;
                crack = crackInTwo(points,*distances,T->start,T->end,epsilon, query, result, ids, &maxDistance);

                T->rec = query;
                T->epsilon = epsilon;
                if(crack >= T->start && T->end >= crack+1)
                {
                    T = T->insertLeft(T, NULL, T->start, crack,queryId);
                    T = T->insertRight(T, NULL, crack+1, T->end,queryId );
                }
            }
            else
            {
                dist = distance(T->rec,query,dimension);
                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++) 
                        {
                            (*result)++;
                        }
                        searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                    }
                    else
                    {
                        searchAndCrack(T->left, query, epsilon, points, distances, result, ids, dimension, queryId);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                        }
                    }
                }
                else
                {
                    searchAndCrack(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                }
            }
        }
    }

    namespace standardMediocre
    {

        int crackInTwoMediocre(float **points, float *distances, int start, int end, float epsilon, float newE, float* query, int *found_count, int *ids, float *maxDistance)
        {
            int i = start, j = end;
            while(true)
            {
                while(distances[i] <= newE && i < end+1){
                    if(distances[i] <= epsilon)
                    {
                        if(distances[i] > *maxDistance)
                            (*maxDistance) = distances[i];
                        (*found_count)++;
                    }
                    i++;
                }

                while(distances[j] > newE && j>start-1){
                    if(distances[j] <= epsilon)
                    {
                        if(distances[j] > *maxDistance)
                            (*maxDistance) = distances[j];
                        (*found_count)++;
                    } 
                    j--;   
                }
                if(i>=j)
                    break;
                swap(points[i], points[j]);       
                swap(distances[i], distances[j]);
                swap(ids[i], ids[j]);
            }

            return j;
        }

        void searchAndCrackMediocre(Node *T, float *query, float epsilon, float **points, float **distances, int *result, int *ids, int dimension, int queryId)
        {   
            int i, j, crack;
            float maxDistance, dist, median;
            vector<int> rnd_dists;

            if(T->right == NULL  && T->left == NULL)
            {
                for(i = T->start; i <= T->end; i++)
                    (*distances)[i] = distance(points[i],query,dimension);

                if((T->end - T->start + 1) <= CRACKTHRESHOLD){
                    for(int i = T->start; i <= T->end ; i++){
                        if((*distances)[i] <= epsilon){
                            (*result)++;
                        }
                    }
                }
                else
                {
                    rnd_dists.clear();
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));

                    median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));

                    maxDistance = -1.0;
                    crack = crackInTwoMediocre(points, *distances, T->start,T->end,epsilon, median, query, result, ids, &maxDistance);

                    T->rec = query;
                    T->epsilon = median;
                    
                    if(crack >= T->start && T->end >= crack+1)
                    {
                        T = T->insertLeft(T, NULL, T->start, crack,queryId);
                        T = T->insertRight(T, NULL, crack+1, T->end,queryId );
                    }
                }
            }
            else
            {
                dist = distance(T->rec,query,dimension);
                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++)
                        {
                            (*result)++;
                        }
                        searchAndCrackMediocre(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                    }
                    else
                    {
                        searchAndCrackMediocre(T->left, query, epsilon, points, distances, result, ids, dimension, queryId);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrackMediocre(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                        }
                    }
                }
                else
                {
                    searchAndCrackMediocre(T->right, query, epsilon, points, distances, result, ids, dimension, queryId);
                }
            }
        }

    }

    namespace caching{

        int partition(int low, int high, int* lp) 
        { 
            if (distances[low] > distances[high]) 
            {
                swap(ids[low], ids[high]); 
                swap(distances[low], distances[high]); 
                swap(points[low], points[high]); 
            }
            int j = low + 1; 
            int g = high - 1, k = low + 1;
            float p = distances[low], q = distances[high]; 
            while (k <= g) { 
        
                if (distances[k] < p) { 
                    swap(ids[k], ids[j]); 
                    swap(distances[k], distances[j]); 
                    swap(points[k], points[j]); 
                    j++; 
                } 
        
                else if (distances[k] >= q) { 
                    while (distances[g] > q && k < g) 
                        g--; 
                    swap(ids[k], ids[g]); 
                    swap(distances[k], distances[g]); 
                    swap(points[k], points[g]); 
                    g--; 
                    if (distances[k] < p) { 
                        swap(ids[k], ids[j]);
                        swap(distances[k], distances[j]);
                        swap(points[k], points[j]);
                        j++; 
                    } 
                } 
                k++; 
            } 
            j--; 
            g++; 
        
            swap(ids[low], ids[j]); 
            swap(ids[high], ids[g]); 

            swap(distances[low], distances[j]); 
            swap(distances[high], distances[g]);

            swap(points[low], points[j]); 
            swap(points[high], points[g]);

            *lp = j;
        
            return g; 
        } 

        void DualPivotQuickSort(int low, int high) 
        { 
            if (low < high) { 
                int lp, rp; 
                rp = partition(low, high, &lp); 
                DualPivotQuickSort(low, lp - 1); 
                DualPivotQuickSort(lp + 1, rp - 1); 
                DualPivotQuickSort(rp + 1, high); 
            } 
        } 

        int binarySearchFirstBigger(float **R, int start, int end, float key)
        {
            if(start >= end)
                return 0;
            
            int low = start;
            int high = end;
            
            int ans = end;
            bool flag = false;
            while (low <= high) {
                
                int mid = (high+low)/2;
                double midVal = (*R)[mid];
                
                if (midVal < key) {
                    low = mid + 1;
                }
                else if (midVal >= key) {
                    ans = mid;
                    high = mid - 1;
                }
            }
            if(low > end)
                return end+1;
            return ans;
        }

        int binarySearchLastSmaller(float **R, int start, int end, float key)
        {
            if(start >= end)
                return 0;
            
            int low = start;
            int high = end;
            
            int ans = end;
            bool flag = false;
            while (low <= high) {
                
                int mid = (high+low)/2;
                double midVal = (*R)[mid];
                
                if (midVal <= key) {
                    low = mid + 1;
                    ans = mid;
                }
                else if (midVal > key) {
                    high = mid - 1;
                }
            }
            if(high < start)
                return start-1;
            return ans;
        }

        void searchAndCrackMediocreCache(Node *T, float *query, float epsilon, float **points, float **distances, int *result, int *ids, int dimension, int queryId, Node* prev, float qDist)
        {   
            int i, j, crack;
            float maxDistance, dist, median;
            vector<int> rnd_dists;
            int low, high;

            
            if(T->right == NULL  && T->left == NULL)
            {
                if((T->end - T->start + 1) <= CRACKTHRESHOLD)
                {
                    auto begin = clock();
                    float qDist = distance(prev->rec, query, dimension);
                    if(T->leftchild == 1)
                    {
                        if(qDist > epsilon)                                                                 // CASE L1
                        {
                            if((*distances)[T->start] >= qDist - epsilon)
                                low = T->start;
                            else
                                low = binarySearchFirstBigger(distances, T->start, T->end, qDist-epsilon);

                            if(qDist >= prev->epsilon)  
                                high = T->end;
                            else
                            {
                                high = binarySearchLastSmaller(distances, T->start, T->end, qDist + epsilon);
                            }

                            for(i = low; i <= high; i++)
                            {
                                if( distance(points[i],query,dimension) <= epsilon)
                                {
                                    (*result)++;
                                }
                            }
                        }
                        else            // CASE L2
                        {
                            if(qDist + prev->epsilon <= epsilon)
                                low = T->end+1;
                            else
                                low = binarySearchFirstBigger(distances, T->start, T->end, epsilon-qDist);

                            for(i = T->start; i < low; i++)
                            {
                                (*result)++;
                            }
                            high = T->end;

                            for(i = low; i <= T->end; i++)
                            {
                                if(distance(points[i],query,dimension) <= epsilon)
                                {
                                    (*result)++;
                                }
                            }
                        }
                    } 
                    else
                    {
                        if(qDist > epsilon+prev->epsilon)
                        {
                            // cout << "CASE R1\n";
                            if((*distances)[T->start] >= qDist + epsilon)
                                low = T->start;
                            else
                                low = binarySearchFirstBigger(distances, T->start, T->end, qDist-epsilon);
                        }
                        else
                        {
                            // cout << "CASE R2\n";
                            low = T->start;
                        }
                        

                        if((*distances)[T->end] <= qDist + epsilon)
                            high = T->end;
                        else
                            high = binarySearchLastSmaller(distances, T->start, T->end, qDist+epsilon);

                        for(i = low; i <= high; i++)
                        {
                            if(distance(points[i],query,dimension) <= epsilon)
                            {
                                (*result)++;
                            }  
                        }
                    }
                }
                else{
                
                    for(i = T->start; i <= T->end; i++)
                        (*distances)[i] = distance(points[i],query,dimension);

                    rnd_dists.clear();
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));
                    rnd_dists.push_back(T->start + ( std::rand() % ( T->end - T->start + 1 ) ));

                    median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));
                    
                    maxDistance = -1.0;
                    crack = standardMediocre::crackInTwoMediocre(points, *distances, T->start,T->end,epsilon, median, query, result, ids, &maxDistance);

                    T->rec = query;
                    T->epsilon = median;


                    if(crack >= T->start && T->end >= crack+1)
                    {
                        T->id = queryId;
                        T = T->insertLeft(T, NULL, T->start, crack,queryId, 1);
                        T = T->insertRight(T, NULL, crack+1, T->end,queryId, 0 );
                        
                        if(crack - T->start + 1 <= CRACKTHRESHOLD)
                        {
                            DualPivotQuickSort(T->start, crack );  
                        }

                        if(T->end - crack <= CRACKTHRESHOLD)
                        {
                            DualPivotQuickSort(crack+1, T->end);  
                        }

                    }
                }
            }
            else
            {
                dist = distance(T->rec,query,dimension);

                if(dist < epsilon + T->epsilon)
                {
                    if(dist < epsilon - T->epsilon)
                    {
                        for(j = T->left->start; j <= T->left->end; j++) 
                        {
                            (*result)++;
                        }
                        searchAndCrackMediocreCache(T->right, query, epsilon, points, distances, result, ids, dimension, queryId, T, dist);
                    }
                    else
                    {
                        searchAndCrackMediocreCache(T->left, query, epsilon, points, distances, result, ids, dimension, queryId, T, dist);
                        if(dist >= T->epsilon - epsilon)
                        {
                            searchAndCrackMediocreCache(T->right, query, epsilon, points, distances, result, ids, dimension, queryId, T, dist);
                        }
                    }
                }
                else
                {
                    searchAndCrackMediocreCache(T->right, query, epsilon, points, distances, result, ids, dimension, queryId, T, dist);
                }
            }
        }    
    }

#endif


int main(int argc, char **argv) {

    char c;
    int runProcessingMethod = -1;
    int res=0;
    float *radius, dist;
    int numPoints=0, dimension=0, numQueries;
    clock_t begin,end, beginPerQuery, endPerQuery;

    #ifdef L2
        std::cout << "L2 distance function" << std::endl;   
    #elif ED
        std::cout << "Edit Distance function" << std::endl;   
    #elif L1
        std::cout << "L1 distance function" << std::endl;   
    #endif

    while ((c = getopt(argc, argv, "lsmcht:")) != -1)
    {
        switch (c)
        {
            case 't':
                CRACKTHRESHOLD = atoi(optarg);
               break;
            case 'l':
                runProcessingMethod = LINEAR_SCAN;
                break;
            case 's':
                runProcessingMethod = STANDARD;
                break;
            case 'm':
                runProcessingMethod = STANDARD_MEDIOCRE;
                break;
            case 'c':
                runProcessingMethod = MEDIOCRE_CACHE;
                break;
            case '?':
                exit(0);
            case 'h':
            default:                
                break;
        }
    }
    cout << "Point dataset : " << argv[optind] << endl;
    cout << "Query dataset : " << argv[optind+1] << endl;
    cout << "threshold: " << CRACKTHRESHOLD << endl;
    #ifdef ED
        words = tools::loadDataString(argv[optind],&numPoints);
        tools::loadQueriesString(argv[optind+1],&numQueries,&radius);
    #else
        tools::loadData(argv[optind], &points, &numPoints, &dimension);
        tools::loadQueries(argv[optind+1], &queries, &radius, &numQueries, dimension );
    #endif 

    distances = (float *)calloc(numPoints,sizeof(float));	
    ids = (int *)malloc(numPoints*sizeof(int));
    for(int i = 0; i < numPoints; i++)
    {
        ids[i] = i;
    }


    #ifdef ED
        switch (runProcessingMethod)
        {
            case LINEAR_SCAN:
                cout << "Linear" << endl;
                cout << "Results\tQueryTime"<<endl;
                begin = clock();
                
                for(int i = 0; i < numQueries; i++)
                {
                    beginPerQuery = clock();
                    unsigned long long result = 0;
                    for( int j = 0; j < numPoints; j++)
                    {
                        // 
                        dist = distance(words[j],(char*)query_words[i].c_str(),0);
                        if(dist <= radius[i])
                        {
                            result++;
                        } 
                    }
                    
                    endPerQuery = clock() - beginPerQuery;
                    cout << result << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;

                break;

            case STANDARD:
                cout << "Standard" << endl;
                cout << "Results\tQueryTime"<<endl;

                root = root->insertString(root, 1.0, 0, numPoints-1,-1 );
                begin = clock();
                for(int i = 0; i < numQueries; i++)
                {
                    res = 0;
                    beginPerQuery = clock();
                    standard::searchAndCrack(root, (char *)query_words[i].c_str(), radius[i], words, &distances, &res, ids, dimension,i);
                    endPerQuery = clock() - beginPerQuery;
                    cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
                break;

            case MEDIOCRE_CACHE:
                cout << "Standard mediocre Cache" << endl;
                cout << "Results\tQueryTime"<<endl;

                root = root->insertString(root, 1.0, 0, numPoints-1,-1 );
                // numberOfNodes++;
                begin = clock();
                for(int i = 0; i < numQueries; i++)
                {
                    res = 0;
                    beginPerQuery = clock();
                    standard_cache_sort::searchAndCrackMediocreThresholdSortBinaryOpt(root, (char*)query_words[i].c_str(), radius[i], words, &distances, &res, ids, dimension,i, root, 0.0);
                    endPerQuery = clock() - beginPerQuery;
                    cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;

                break;

            case STANDARD_MEDIOCRE:
                cout << "Standard Mediocre Threshold" << endl;
                cout << "Results\tQueryTime"<<endl;

                root = root->insertString(root, 1.0, 0, numPoints-1,-1 );
                // numberOfNodes++;
                begin = clock();
                for(int i = 0; i < numQueries; i++)
                {
                    res = 0;
                    beginPerQuery = clock();
                    standardMediocre::searchAndCrackMediocreThreshold(root, (char*)query_words[i].c_str(), radius[i], words, &distances, &res, ids, dimension,i);
                    endPerQuery = clock() - beginPerQuery;
                    cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;

                break; 

        }
   
    #else
    
    switch (runProcessingMethod)
    {
        case LINEAR_SCAN:
            cout << "Linear" << endl;
            cout << "Results\tQueryTime"<<endl;
            begin = clock();
            
            for(int i = 0; i < numQueries; i++)
            {
                beginPerQuery = clock();
                unsigned long long result = 0;
                for( int j = 0; j < numPoints; j++)
                {
                    dist = distance(points[j],queries[i],dimension);
                    if(dist <= radius[i])
                    {
                        result++;
                    } 
                }
                endPerQuery = clock() - beginPerQuery;
                cout << result << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
            }
            end = clock() - begin;
            cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
            break;

        case STANDARD:
            cout << "Standard" << endl;
            cout << "Results\tQueryTime"<<endl;

            root = root->insert(root, 1.0, 0, numPoints-1,-1 );
            begin = clock();
            for(int i = 0; i < numQueries; i++)
            {
                res = 0;
                beginPerQuery = clock();
                standard::searchAndCrack(root, queries[i], radius[i], points, &distances, &res, ids, dimension,i);
                endPerQuery = clock() - beginPerQuery;
                cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
            }
            end = clock() - begin;
            cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
            break; 

        case STANDARD_MEDIOCRE:
            cout << "Standard Mediocre" << endl;
            cout << "Results\tQueryTime"<<endl;

            root = root->insert(root, 1.0, 0, numPoints-1,-1 );
            begin = clock();
            for(int i = 0; i < numQueries; i++)
            {
                res = 0;
                beginPerQuery = clock();
                standardMediocre::searchAndCrackMediocre(root, queries[i], radius[i], points, &distances, &res, ids, dimension,i);
                endPerQuery = clock() - beginPerQuery;
                cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
            }
            end = clock() - begin;
            cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
            break; 

        case MEDIOCRE_CACHE:
            cout << "Standard mediocre Cache" << endl;
            cout << "Results\tQueryTime"<<endl;

            root = root->insert(root, 1.0, 0, numPoints-1,-1 );
            begin = clock();
            for(int i = 0; i < numQueries; i++)
            {
                res = 0;
                beginPerQuery = clock();
                caching::searchAndCrackMediocreCache(root, queries[i], radius[i], points, &distances, &res, ids, dimension,i, root, 0.0);
                endPerQuery = clock() - beginPerQuery;
                cout << res << "\t" << (double)endPerQuery/CLOCKS_PER_SEC << endl;
            }
            end = clock() - begin;
            cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
            break;
    }

    #endif
}

