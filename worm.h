#ifndef WORM_H
#define WORM_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct book_t {
	size_t id;
	size_t author_id;
	size_t publisher_id;
	size_t* b_author_edges;
	size_t* b_citation_edges;
	size_t* b_publisher_edges;
	size_t n_author_edges;
	size_t n_citation_edges;
	size_t n_publisher_edges;
} book_t;

typedef struct result_t {
	book_t** elements;
	size_t n_elements;
} result_t;

typedef struct bfs_node {
	book_t * book;
	size_t level;
	struct bfs_node * parent;
	size_t index;
} bfs_node;

typedef struct bfs_queue {
	bfs_node** nodes;
	size_t size;
	size_t cap;
	size_t begin;
} bfs_queue;

typedef struct indx{
	size_t index;
	bool seen;
}indx;


typedef struct thread_data {
	int thread_id;
	struct book_t* nodes;
	size_t count;
	size_t id;
	// 1 search by id
	// 2 search by author
	// 3 search by publisher
	size_t mode;
	uint16_t distance;
	book_t * target_book;
	result_t *result;
	bfs_node * bfs_node;
	bool is_seen;
} TDATA;

void* find_worker(void* args);

result_t* find_book(book_t* nodes, size_t count, size_t book_id);
result_t* find_books_by_author(book_t* nodes, size_t count, size_t author_id);
result_t* find_books_reprinted(book_t* nodes, size_t count, size_t publisher_id);
result_t* find_books_k_distance(book_t* nodes, size_t count, size_t book_id, uint16_t k);
result_t* find_shortest_distance(book_t* nodes, size_t count, size_t b1_id, size_t b2_id);
result_t* find_shortest_edge_type(book_t* nodes, size_t count, size_t a1_id, size_t a2_id);

#endif

