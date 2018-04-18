#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include "worm.h"

//global variables
size_t g_nthreads = 4;
bfs_queue *k_queue;
bfs_queue *k_bfs_nodes;



//initialise result
result_t * new_result(){
    result_t * result = malloc(sizeof(result_t));
    result->elements = malloc(sizeof(book_t *));
    result->n_elements = 0;
    return result;
}

//free results after merged
void free_result(result_t * result){
	if (result == NULL){
        return;
    }
    if (result->elements != NULL){
        free(result->elements);
    }
    free(result);
}

int getStart(TDATA *tdatap){
    int items = tdatap->count / g_nthreads; 
    int start = tdatap->thread_id * items;
    return start;
}

int getEnd(TDATA *tdatap){
    int items = tdatap->count / g_nthreads; 
    int end = (tdatap->thread_id + 1) * items;
    if (tdatap->thread_id == g_nthreads - 1) {
        end = tdatap->count;
    }
    return end;
}

//append book to results
void result_append(result_t * result, book_t* node){
    result->elements = realloc(result->elements, sizeof(book_t *) * (result->n_elements + 1));
    result->elements[result->n_elements++] = node;
}

//result merging
void result_append_all(result_t * merged_result, result_t * result){
    if (result->n_elements != 0){
        // merge thread result into result
        merged_result->elements = realloc(merged_result->elements, sizeof(book_t *) * (merged_result->n_elements + result->n_elements));
        
        // copy all the results over
        for (int i = 0; i < result->n_elements; i++){
            merged_result->elements[merged_result->n_elements++] = result->elements[i];
        }
    }
    free_result(result);
}

//worker function searched for matching id depending which mode
void* find_worker(void* args){
    TDATA *tdatap = (TDATA *)args;

    for (int i = getStart(tdatap); i < getEnd(tdatap); i++){
        book_t *node = tdatap->nodes+i;
        bool isMatch = false;

        // compare book id
        if (tdatap->mode == 1 && node->id == tdatap->id) {
            isMatch = true;
        }
        // compare author id
        if (tdatap->mode == 2 && node->author_id == tdatap->id) {
            isMatch = true;
        }
        // compare publisher id
        if (tdatap->mode == 3 && node->publisher_id == tdatap->id) {
            isMatch = true;
        }
        //adds to result
        if (isMatch) {
            result_append(tdatap->result, node);
        }
    }
    return tdatap;
}
//creates threads for finding book by any id
struct result_t* find_by_any_id(struct book_t* nodes, size_t count, size_t id, int mode) {
    
    pthread_t th[g_nthreads];
    TDATA *tdata[g_nthreads];
    //creates thread and calls worker function
    for (int i = 0; i < g_nthreads; ++i) 
    {
        tdata[i] = (TDATA*)malloc(sizeof(TDATA));
        TDATA *tdatap = tdata[i];
        tdatap->thread_id = i;
        tdatap->nodes = nodes;
        tdatap->count = count;
        tdatap->id = id;
        tdatap->mode = mode;
        tdatap->result = new_result();
        pthread_create(th+i, NULL, find_worker, tdatap);
    }
    
    result_t* result = new_result();
    //join threads and merge results
    for (int i = 0; i < g_nthreads; ++i) {
        TDATA *tdatap;
        pthread_join(th[i], (void**)&tdatap);
        result_append_all(result, tdatap->result);
        free(tdatap);
    }
    return result;
}

// Returns result set containing book with given book id.
result_t* find_book(book_t* nodes, size_t count, size_t book_id) {
    result_t * merged_result = find_by_any_id(nodes, count, book_id, 1);
    result_t *result = new_result();
    if(merged_result->n_elements != 0) {
        result_append(result, merged_result->elements[0]);
    }

    free_result(merged_result); 
    return result;
}

// Returns result set containing books by given author id.
result_t* find_books_by_author(book_t* nodes, size_t count, size_t author_id) {
    return find_by_any_id(nodes, count, author_id, 2);
}


// Returns result set containing books that have been reprinted by a different publisher.
result_t* find_books_reprinted(book_t* nodes, size_t count, size_t publisher_id) {
    // all the books belongs to this publisher
    result_t * result_books = find_by_any_id(nodes, count, publisher_id, 3);
	
    result_t * final_result = new_result();
    int n_ele = result_books->n_elements;
	for(int i = 0; i < count; i++){
		for(int j = 0; j < n_ele; j++){
			if (nodes[i].id == result_books->elements[j]->id
            && nodes[i].author_id == result_books->elements[j]->author_id
            && nodes[i].publisher_id != result_books->elements[j]->publisher_id) {
            	result_append(final_result, &nodes[i]);
       		}	
		}
	}
	
	
    free_result(result_books);
    
    return final_result;
}

//find books with same id for k distance, 
book_t * find_book_by_id(book_t* nodes, size_t count, size_t book_id, int* index_ptr) {
    for (int i = 0; i < count; i++) {
        if (nodes[i].id == book_id) {
            *index_ptr = i;
        }
    }
    if (*index_ptr == -1) {
        return NULL;
    } else {
        return &nodes[*index_ptr];
    }
}

result_t * find_books_by_id(book_t* nodes, size_t count, size_t book_id) {
    result_t * books = new_result();
    for (int i = 0; i < count; i++) {
        if (nodes[i].id == book_id) {
            result_append(books, &nodes[i]);
        }
    }
    return books;
}

//adds to the queue
void enqueue(bfs_queue *queue, bfs_node* node) {
    //first time enqueue
    if(queue->nodes == NULL){
        queue->cap = queue->cap *2;
        queue->nodes = malloc(sizeof(bfs_node*)*queue->cap);
    } else if(queue->size == queue->cap){
        queue->cap = queue->cap *2;
        queue->nodes = realloc(queue->nodes, sizeof(bfs_node*) * queue->cap);
    }
    queue->nodes[queue->size] = node;
    queue->size++;
}

//removes fron the front of queue
bfs_node* dequeue(bfs_queue *queue) {
    if (queue->size - queue->begin == 0) {
        return NULL;
    }
    bfs_node* node = queue->nodes[queue->begin];
    queue->begin++;
    return node;
}

//creates bfs nodes
bfs_node * new_bfs_node(book_t * book, size_t level, size_t index) {
    bfs_node *root = malloc(sizeof(bfs_node));
    root->book = book;
    root->level = level;
    root->parent = NULL;
    root->index = index;
    return root;
}

//nodes with parent to keep trace back the path
bfs_node * new_bfs_node_with_parent(book_t * book, size_t level, bfs_node* parent, size_t index) {
    bfs_node *root = malloc(sizeof(bfs_node));
    root->book = book;
    root->level = level;
    root->parent = parent;
    root->index = index;
    return root;
}

//worker function for is_book_seen
void* is_book_seen_worker(void* args){
    TDATA *tdatap = (TDATA *)args;

    int items = tdatap->result->n_elements / g_nthreads; 
    int start = tdatap->thread_id * items;
    int end = (tdatap->thread_id + 1) * items;
    if (tdatap->thread_id == g_nthreads - 1) {
        end = tdatap->result->n_elements;
    }
    for (int i = start; i < end; i++) {
        if (tdatap->result->elements[i] == tdatap->target_book ) {
            tdatap->is_seen = true;
            break;
        }
    }

    return tdatap;
}

//checks if book has appeared in queue, creates thread
bool is_book_seen(result_t * result, book_t* book) {
	
    //implementation without multi-thread
    for (int i = 0; i < result->n_elements; i++) {
        if (result->elements[i]->id == book->id) {
            return true;
        }
    }
    return false;
}

//function to create new bfs queue
bfs_queue * new_queue() {
    
    bfs_queue *queue = malloc(sizeof(bfs_queue));
    //initialise size and cap
    queue->size = 0;
    queue->cap = 100;
    queue->begin = 0;
    queue->nodes = NULL;
    return queue;
}


// Returns result set containing books that are k distance from given book.
result_t* find_books_k_distance(book_t* nodes, size_t count, size_t book_id, uint16_t k) {

    // keep the seen one
    result_t * result = new_result();
    int found_i[count];
    for (int i = 0; i < count; i++) {
        found_i[i] = 0;
    }

    k_queue = new_queue();
    k_bfs_nodes = new_queue();


    int index1 = -1;
    book_t *book = find_book_by_id(nodes, count, book_id, &index1);
    if (book == NULL) {
        return result;
    }

  
    bfs_node *root = new_bfs_node(book, 0, index1);
    enqueue(k_queue, root);
    enqueue(k_bfs_nodes, root);

    bool found = false;

    while(k_queue->size - k_queue->begin > 0) {

        bfs_node * u = dequeue(k_queue);
        
        if (u->level > k) {
            found = true;
            break;
        }

        if (found_i[u->index] == 0) {
            found_i[u->index] = 1;
            result_append(result, u->book);
          
            //implementation without multi-thread
            for (int i = 0; i < u->book->n_citation_edges; i++) {
                size_t citation_index = u->book->b_citation_edges[i];
                if (citation_index < count) {
                    book_t *neighbour_book = &nodes[citation_index];
                    bfs_node *w = new_bfs_node(neighbour_book, u->level + 1, citation_index);
                    enqueue(k_queue, w);
                    enqueue(k_bfs_nodes, w);
                }
            }

        }

    }

    // free the rest elements of queue
    if (k_bfs_nodes->nodes != NULL) {
        for (int i = 0; i < k_bfs_nodes->size; i++) {
            free(k_bfs_nodes->nodes[i]);
        }
        free(k_bfs_nodes->nodes);
    }
    free(k_bfs_nodes);

    if (k_queue->nodes != NULL) {
        free(k_queue->nodes);
    }
    free(k_queue);

    //reset results when nothing is found. return empty results
    return result;
   
}


// Returns result set containing books in shortest path between book 1 and 2.
result_t* find_shortest_distance(book_t* nodes, size_t count, size_t b1_id, size_t b2_id) {

    int index1 = -1;
    book_t *book = find_book_by_id(nodes, count, b1_id, &index1);
    int found_i[count];
    for (int i = 0; i < count; i++) {
        found_i[i] = 0;
    }

    if (book == NULL) {
        return new_result();
    }
    int index2 = -1;
    book_t *book2 = find_book_by_id(nodes, count, b2_id, &index2);
    if (book2 == NULL) {
        return new_result();
    }


    bfs_node *root = new_bfs_node_with_parent(book, 0, NULL, index1);
    k_queue = new_queue();
    k_bfs_nodes = new_queue();


    
    enqueue(k_queue, root);
    enqueue(k_bfs_nodes, root);
    bfs_node *dest = NULL;

    // mark other nodes with same book_id as seen as well
    for (int i = 0; i < count; i++) {
        if (nodes[i].id == b1_id && &nodes[i] != book) {
            found_i[i]++;
        }
    }


    bool found = false;

    while(k_queue->size != k_queue->begin) {
        bfs_node * u = dequeue(k_queue);
        
        if (found) {
            break;
        }

        //treat all edges as one tree
        
        if (!found_i[u->index]) {
            found_i[u->index]++;
            
            for (int i = 0; i < u->book->n_citation_edges; i++) {
                size_t index = u->book->b_citation_edges[i];
                book_t * neighbour_book = &nodes[index];

                if (neighbour_book != NULL) {
                    bfs_node *w = new_bfs_node_with_parent(neighbour_book, u->level + 1, u, index);
                    enqueue(k_queue, w);
                    enqueue(k_bfs_nodes, w);
                    if (w->book->id == b2_id) {
                        found = true;
                        dest = w;
                        break;
                    }
                }
            }
			
            //loops through publisher edges using bfs 
            for (int i = 0; i < u->book->n_publisher_edges; i++) {
                size_t index = u->book->b_publisher_edges[i];
                book_t * neighbour_book = &nodes[index];

                if (neighbour_book != NULL) {
                    bfs_node *w = new_bfs_node_with_parent(neighbour_book, u->level + 1, u, index);
                    enqueue(k_queue, w);
                    enqueue(k_bfs_nodes, w);
                    if (w->book->id == b2_id) {
                        found = true;
                        dest = w;
                        break;
                    }
                }
            }
            //loops through author edges using bfs 
            for (int i = 0; i < u->book->n_author_edges; i++) {
                size_t index = u->book->b_author_edges[i];
                book_t * neighbour_book = &nodes[index];

                if (neighbour_book != NULL) {
                    bfs_node *w = new_bfs_node_with_parent(neighbour_book, u->level + 1, u, index);
                    enqueue(k_queue, w);
                    enqueue(k_bfs_nodes, w);
                    if (w->book->id == b2_id) {
                        found = true;
                        dest = w;
                        break;
                    }
                }
            }

        }

    }

    //trace back to get path 
    result_t * path = new_result();
    bfs_node *n = dest;
    while (n != NULL) {
        result_append(path, n->book);
        n = n->parent;
    }
    

    // free the rest elements of k_queue
    if (k_bfs_nodes->nodes != NULL) {
        for (int i = 0; i < k_bfs_nodes->size; i++) {
            free(k_bfs_nodes->nodes[i]);
        }
        free(k_bfs_nodes->nodes);
    }
    free(k_bfs_nodes);

    if (k_queue->nodes != NULL) {
        free(k_queue->nodes);
    }
    free(k_queue);

    if(path->n_elements == 0){
        free_result(path);
        result_t * path2 = new_result();
        return path2;
    }
    
  
    return path;
}






