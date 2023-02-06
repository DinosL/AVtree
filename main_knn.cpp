#include "def.h"

int CRACKTHRESHOLD=0, K = 0,  biggestWordSize = 0;
Node *root;
vector<string> query_words;
char **words;
int *ids, *gc;

#ifdef L2
    #define distance(point,query,dimension) tools::L2Distance(point,query,dimension)
    typedef tuple<float, Node *, int> knn_ext_pair;
    typedef pair<float, int> knn_pair;
#elif L1
     #define distance(point,query,dimension) tools::L1Distance(point,query,dimension)
     typedef tuple<float, Node *, int> knn_ext_pair;
    typedef pair<float, int> knn_pair;
#elif ED
    #define distance(point,query,dimension) (tools::EditDistance(point,query))
    typedef tuple<float, Node *, char*> knn_ext_pair;   
    typedef pair<float, char *> knn_pair;
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

    void showpq(priority_queue<knn_ext_pair> gq)
    {
        priority_queue<knn_ext_pair> g = gq;
        while (!g.empty()) {
            cerr << get<2>(g.top()) << " "; 
            g.pop();
        }
        cerr << '\n';
    }

    void showpq(priority_queue<knn_pair> gq)
    {
        priority_queue<knn_pair> g = gq;
        while (!g.empty()) {
            cerr << g.top().second << " ";
            g.pop();
        }
        cerr << '\n';
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
            // char *ptr1 = strtok(line, ",");
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
    namespace TwoPQStrings {
        int crackInTwo(char **points, float *distances, int start, int end, float epsilon, int *ids)
        {
            int i = start, j = end;
            while(true)
            {
                while(distances[i] <= epsilon && i < end+1){
                    i++;
                }
                while(distances[j] > epsilon && j>start-1){
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

        Node* knnMediocre(Node *T,char *query, char **points, float **distances, priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> &guide_pq, priority_queue<knn_ext_pair> &result, int *ids, int dimension, int queryId)
        {
            int actualK = K;
            int i, ii, crack;
            float median, leftMinDist, rightMinDist, dist, distTop;
            vector<int> rnd_dists;
            Node *topElement;

            while(!guide_pq.empty() && ( result.size() < K || get<0>(guide_pq.top()) < get<0>(result.top())))
            {
                Node* node = get<1>(guide_pq.top());
                if(node->right == NULL && node->left == NULL)
                {
                    for(i = node->start; i <= node->end; i++)
                    {
                        
                        (*distances)[i] = distance(points[i],query,0);
                        if(result.size() < K)
                        {
                            result.push(make_tuple((*distances)[i],node,words[i]));
                        }
                        else if( (*distances)[i] < get<0>(result.top()))
                        {
                            result.pop();
                            result.push(make_tuple((*distances)[i],node,words[i]));
                        }
                        
                    }
                    guide_pq.pop();
                    if((node->end - node->start + 1) > CRACKTHRESHOLD)
                    {
                        rnd_dists.clear();
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));

                        median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));
                        crack = crackInTwo(points, *distances, node->start, node->end, median, ids);

                        node->rec = query;
                        node->epsilon = median;
                        
                        if(crack >= node->start && node->end >= crack+1)
                        {
                            node = node->insertLeftString(node, NULL, node->start, crack,queryId);
                            node = node->insertRightString(node, NULL, crack+1, node->end,queryId );
                        }
                    }
                }
                else
                {
                    dist = distance(node->rec, query,0);
                    guide_pq.pop();

                    leftMinDist = max(0.0f,dist-node->epsilon);
                    rightMinDist =  max(0.0f,node->epsilon-dist);
                    if(result.size() < K)
                    {
                        guide_pq.push(make_tuple(leftMinDist,node->left,words[0]));
                        guide_pq.push(make_tuple(rightMinDist,node->right,words[0]));
                    }
                    else{
                        if(leftMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(leftMinDist,node->left,words[0]));
                        }
                        
                        if(rightMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(rightMinDist,node->right,words[0]));
                        }
                    }
                }
            }
            return T;

        }
}



#else
    namespace TwoPQ {
        int crackInTwo(float **points, float *distances, int start, int end, float epsilon, int *ids)
        {
            int i = start, j = end;
            while(true)
            {
                while(distances[i] <= epsilon && i < end+1){
                    i++;
                }
                while(distances[j] > epsilon && j>start-1){
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

        Node* knnMediocre(Node *T,float *query, float **points, float **distances, priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> &guide_pq, priority_queue<knn_ext_pair> &result, int *ids, int dimension, int queryId)
        {
            int actualK = K;
            int i, ii, crack;
            float median, leftMinDist, rightMinDist, dist, distTop;
            vector<int> rnd_dists;
            Node *topElement;

            while(!guide_pq.empty() && ( result.size() < K || guide_pq.top() < result.top()))
            {
                Node* node = get<1>(guide_pq.top());
                if(node->right == NULL && node->left == NULL)
                {
                    for(i = node->start; i <= node->end; i++)
                    {
                        (*distances)[i] = distance(points[i],query,dimension);
                        if(result.size() < K)
                        {
                            result.push(make_tuple((*distances)[i],node,ids[i]));
                        }
                        else if( (*distances)[i] < get<0>(result.top()))
                        {
                            result.pop();
                            result.push(make_tuple((*distances)[i],node,ids[i]));
                        }
                    }
                    guide_pq.pop();
                    if((node->end - node->start + 1) > CRACKTHRESHOLD)
                    {
                        rnd_dists.clear();
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));
                        rnd_dists.push_back(node->start + ( std::rand() % ( node->end - node->start + 1 ) ));

                        median = max(min((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]), min(max((*distances)[rnd_dists[0]],(*distances)[rnd_dists[1]]),(*distances)[rnd_dists[2]]));
                        crack = crackInTwo(points, *distances, node->start, node->end, median, ids);

                        node->rec = query;
                        node->epsilon = median;
                        
                        if(crack >= node->start && node->end >= crack+1)
                        {
                            node = node->insertLeft(node, NULL, node->start, crack,queryId);
                            node = node->insertRight(node, NULL, crack+1, node->end,queryId );

                        }
                    }
                }
                else
                {
                    dist = distance(node->rec, query, dimension);
                    guide_pq.pop();

                    leftMinDist = max(0.0f,dist-node->epsilon);
                    rightMinDist =  max(0.0f,node->epsilon-dist);
                    if(result.size() < K)
                    {
                        guide_pq.push(make_tuple(leftMinDist,node->left,-1));
                        guide_pq.push(make_tuple(rightMinDist,node->right,-1));
                    }
                    else{
                        if(leftMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(leftMinDist,node->left,-1));
                        }
                        
                        if(rightMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(rightMinDist,node->right,-1));
                        }
                    }
                }
            }
            return T;

        }

        Node* knnCrackOnKth(Node *T,float *query, float **points, float **distances, priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> &guide_pq, priority_queue<knn_ext_pair> &result, int *ids, int dimension, int queryId)
        {
            int actualK = K, ii,i;
            float dist, leftMinDist, rightMinDist;
            Node* node;

            while(!guide_pq.empty() && ( result.size() < K || guide_pq.top() < result.top()))
            {
                node = get<1>(guide_pq.top());
                if(node->right == NULL && node->left == NULL)
                {
                    for(i = node->start; i <= node->end; i++)
                            (*distances)[i] = distance(points[i],query,dimension);
                    
                    guide_pq.pop();
                    for(ii = node->start; ii <= node->end; ii++)
                    {
                        if(result.size() < K)
                        {
                            result.push(make_tuple((*distances)[ii],node,ids[ii]));
                        }
                        else if( (*distances)[ii] < get<0>(result.top()))
                        {
                            result.pop();
                            result.push(make_tuple((*distances)[ii],node,ids[ii]));
                        }
                    }
                }
                else
                {
                    dist = distance(node->rec, query, dimension);
                    guide_pq.pop();

                    leftMinDist = max(0.0f,dist-node->epsilon);
                    rightMinDist =  max(0.0f,node->epsilon-dist);

                    if(result.size() < K)
                    {
                        guide_pq.push(make_tuple(leftMinDist,node->left,-1));
                        guide_pq.push(make_tuple(rightMinDist,node->right,-1));
                    }
                    else{
                        if(leftMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(leftMinDist,node->left,-1));
                        }
                        
                        if(rightMinDist < get<0>(result.top()))
                        {
                            guide_pq.push(make_tuple(rightMinDist,node->right,-1));
                        }
                    }
                }
            }

            float kthDist = get<0>(result.top());
            while(!result.empty())
            {
                Node * node = get<1>(result.top());
                cerr << get<2>(result.top()) << " ";
                result.pop();
                if((node->end - node->start + 1) > CRACKTHRESHOLD)
                {
                    int crack = crackInTwo(points, *distances, node->start, node->end, kthDist, ids);

                    node->rec = query;
                    node->epsilon = kthDist;
                    
                    if(crack >= node->start && node->end >= crack+1)
                    {
                        node = node->insertLeft(node, NULL, node->start, crack,queryId);
                        node = node->insertRight(node, NULL, crack+1, node->end,queryId );
                    }
                }
            }
            return T;
        }
    }

#endif


int main(int argc, char **argv) {

    #ifdef L2
        std::cout << "L2" << std::endl;   
    #elif ED
        std::cout << "Edit Distance" << std::endl;   
    #elif L1
        std::cout << "L1" << std::endl;   
    #endif

    char c;
    int runProcessingMethod = -1;
    int res=0,i,j;
    float **points, **queries, *radius, *distances;
    
    int numPoints=0, dimension=0, numQueries;
    clock_t begin,end, beginPerQuery, endPerQuery;

    while ((c = getopt(argc, argv, "lmoht:n:")) != -1)
    {
        switch (c)
        {
            case 't':
                CRACKTHRESHOLD = atoi(optarg);
                break;
            case 'l':
                runProcessingMethod = LINEAR_SCAN;
                break;
            case 'n':
                K = atoi(optarg);
                break;
            case 'm':
                runProcessingMethod = KNNMediocre;
                break;
            case 'o':
                runProcessingMethod = KNNCrackOnKth;
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
    cout << "Find " << K << "nn" << endl;
    if(K==0)
    {
        cout << "nn value in kNN cannot be zero\n";
        exit(1);
    }

    #ifdef ED
        words = tools::loadDataString(argv[optind],&numPoints);
        tools::loadQueriesString(argv[optind+1],&numQueries,&radius);
    #else
        tools::loadData(argv[optind], &points, &numPoints, &dimension);
        tools::loadQueries(argv[optind+1], &queries, &radius, &numQueries, dimension );
    #endif  

    distances = (float *)malloc(numPoints*sizeof(float));	
    ids = (int *)malloc(numPoints*sizeof(int));
    for(int i = 0; i < numPoints; i++)
    {
        ids[i] = i;
    }
    float dist = 0.0;

    #ifdef ED
        switch (runProcessingMethod)
        {
            case LINEAR_SCAN:
        
            cout << "Linear KNN" << endl;
            cout << "QueryTime\tDistanceComputation"<<endl;
            begin = clock();
            for(i = 0; i < numQueries; i++)
            {
                beginPerQuery = clock();
                priority_queue<knn_pair> result;
                for(j = 0; j < numPoints; j++)
                {
                    dist = distance(words[j],(char*)query_words[i].c_str(),0);
                    if(result.size() < K)
                    {
                        result.push(make_pair(dist,words[j]));
                    }
                    else if( dist < result.top().first)
                    {
                        result.pop();
                        result.push(make_pair(dist,words[j]));
                    }
                }
                
                endPerQuery = clock() - beginPerQuery;
                tools::showpq(result);
                // cerr << "\n";
                cout << (double)endPerQuery/CLOCKS_PER_SEC << endl;
            }
            end = clock() - begin;
            cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
            break;

            case KNNMediocre:
                cout << "KNN guide is minPQ and result is maxPQ" << endl;
                cout << "QueryTime" << endl;

                root = root->insertString(root, -1.0, 0, numPoints-1,-1 );

                begin = clock();
                for(int i = 0; i < numQueries; i++)
                {
                    priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> guide_pq;
                    priority_queue<knn_ext_pair> result;

                    guide_pq.push(make_tuple(10,root,words[0]));
                    beginPerQuery = clock();
                    TwoPQStrings::knnMediocre(root,(char*)query_words[i].c_str(),words,&distances,guide_pq,result,ids,dimension,i);
                    endPerQuery = clock() - beginPerQuery;
                    tools::showpq(result);
                    cout << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
                
                break;
        }
    #else
        switch (runProcessingMethod)
        {
            case LINEAR_SCAN:
            
                cout << "Linear KNN" << endl;
                cout << "QueryTime\tDistanceComputation"<<endl;
                begin = clock();
                for(i = 0; i < numQueries; i++)
                {
                    beginPerQuery = clock();
                    priority_queue<knn_pair> result;

                    for(j = 0; j < numPoints; j++)
                    {
                        dist = distance(points[j],queries[i],dimension);
                        if(result.size() < K)
                        {
                            result.push(make_pair(dist,j));
                        }
                        else if( dist < get<0>(result.top()))
                        {
                            result.pop();
                            result.push(make_pair(dist,j));
                        }
                    }
                    
                    endPerQuery = clock() - beginPerQuery;
                    int current_size = 0;

                    while (!result.empty()) {
                        if(current_size == K)
                            break;
                        cerr << result.top().second <<" ";
                        result.pop();
                        current_size++;
                    }
                    cerr << "\n";
                    endPerQuery = clock() - beginPerQuery;
                    cout << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
                break;

            case KNNMediocre:
                cout << "KNN guide is minPQ and result is maxPQ" << endl;
                cout << "QueryTime" << endl;

                root = root->insert(root, -1.0, 0, numPoints-1,-1 );

                begin = clock();
                for(int i = 0; i < numQueries; i++)
                {
                    priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> guide_pq;
                    priority_queue<knn_ext_pair> result;

                    guide_pq.push(make_tuple(10,root,-1));
                    beginPerQuery = clock();
                    TwoPQ::knnMediocre(root,queries[i],points,&distances,guide_pq,result,ids,dimension,i);
                    endPerQuery = clock() - beginPerQuery;
                    tools::showpq(result);
                    cout << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
                
                break;

            case KNNCrackOnKth:
                cout << "KNN crack on K-th" << endl;
                cout << "QueryTime" << endl;

                root = root->insert(root, -1.0, 0, numPoints-1,-1 );

                begin = clock();
                for(i = 0; i < numQueries; i++)
                {
                    priority_queue<knn_ext_pair,vector<knn_ext_pair>,greater<knn_ext_pair>> guide_pq;
                    priority_queue<knn_ext_pair> result;
                    guide_pq.push(make_tuple(10,root,-1));
                    beginPerQuery = clock();
                    TwoPQ::knnCrackOnKth(root,queries[i],points,&distances,guide_pq,result,ids,dimension,i);
                    endPerQuery = clock() - beginPerQuery;
                    tools::showpq(result);
                    cout << (double)endPerQuery/CLOCKS_PER_SEC << endl;
                }
                end = clock() - begin;
                cout << "Total time " << (double)end/CLOCKS_PER_SEC << endl;
                break;
        }

    #endif
}