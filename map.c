#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bseq.h"
#include "kvec.h"
#include "minimap.h"
#include "sdust.h"

void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->radius = 500;
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->sdust_thres = 25;
	opt->flag = 0;
}

/****************************
 * Find approxiate mappings *
 ****************************/

struct mm_tbuf_s { // per-thread buffer
	mm128_v mini; // query minimizers
	mm128_v coef; // Hough transform coefficient
	mm128_v intv; // intervals on sorted coef
	uint32_v reg2mini;
	uint32_v rep_aux;
	sdust_buf_t *sdb;
	// the following are for computing LIS
	uint32_t n, m;
	uint64_t *a;
	size_t *b, *p;
	// final output
	kvec_t(mm_reg1_t) reg;
};

mm_tbuf_t *mm_tbuf_init()
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	b->sdb = sdust_buf_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	free(b->mini.a); free(b->coef.a); free(b->intv.a); free(b->reg.a); free(b->reg2mini.a); free(b->rep_aux.a);
	free(b->a); free(b->b); free(b->p);
	sdust_buf_destroy(b->sdb);
	free(b);
}

#include "ksort.h"
#define sort_key_64(a) (a)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8) 
#define lt_low32(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(low32lt, uint64_t, lt_low32)
#define gt_low32(a, b) ((uint32_t)(a) > (uint32_t)(b))
KSORT_INIT(low32gt, uint64_t, gt_low32)

static void drop_rep(mm_tbuf_t *b, int min_cnt)
{
	int i, j, n, m;
	uint32_t thres;
	b->rep_aux.n = 0;
	for (i = 0; i < b->mini.n; ++i)
		if (b->mini.a[i].y>>32)
			kv_push(uint32_t, b->rep_aux, b->mini.a[i].y>>32);
	if (b->rep_aux.n < 3) return;
	thres = (uint32_t)(ks_ksmall_uint32_t(b->rep_aux.n, b->rep_aux.a, b->rep_aux.n>>1) * MM_DEREP_Q50 + .499);
	for (i = n = m = 0; i < b->reg.n; ++i) {
		int cnt = 0, all_cnt = b->reg.a[i].cnt;
		for (j = 0; j < all_cnt; ++j)
			if (b->mini.a[b->reg2mini.a[m + j]].y>>32 <= thres)
				++cnt;
		if (cnt >= min_cnt)
			b->reg.a[n++] = b->reg.a[i];
		m += all_cnt;
	}
//	printf("%ld=>%d\t%d\n", b->reg.n, n, thres);
	b->reg.n = n;
}

static void proc_intv(mm_tbuf_t *b, int which, int k, int min_cnt, int max_gap)
{
	int i, j, l_lis, rid = -1, rev = 0, start = b->intv.a[which].y, end = start + b->intv.a[which].x;
	if (end - start > b->m) { // make room for arrays needed by LIS
		b->m = end - start;
		kv_roundup32(b->m);
		b->a = (uint64_t*)realloc(b->a, b->m * 8);
		b->b = (size_t*)realloc(b->b, b->m * sizeof(size_t));
		b->p = (size_t*)realloc(b->p, b->m * sizeof(size_t));
	}
	b->n = 0;
	for (i = start; i < end; ++i) // prepare the input array _a_ for LIS
		if (b->coef.a[i].x != UINT64_MAX)
			b->a[b->n++] = b->coef.a[i].y, rid = b->coef.a[i].x << 1 >> 33, rev = b->coef.a[i].x >> 63;
	if (b->n < min_cnt) return;
	radix_sort_64(b->a, b->a + b->n);
	l_lis = rev? ks_lis_low32gt(b->n, b->a, b->b, b->p) : ks_lis_low32lt(b->n, b->a, b->b, b->p); // LIS
	if (l_lis < min_cnt) return;
	for (i = 1, j = 1; i < l_lis; ++i) // squeeze out reused minimizers
		if (b->a[b->b[i]]>>32 != b->a[b->b[i-1]]>>32)
			b->a[b->b[j++]] = b->a[b->b[i]];
	l_lis = j;
	if (l_lis < min_cnt) return;
	for (i = 1, start = 0; i <= l_lis; ++i) { // convert LIS to a region; possibly break LIS at long gaps
		int32_t qgap = i == l_lis? 0 : ((uint32_t)b->mini.a[b->a[b->b[i]]>>32].y>>1) - ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1);
		if (i == l_lis || (qgap > max_gap && abs((int32_t)b->a[b->b[i]] - (int32_t)b->a[b->b[i-1]]) > max_gap)) {
			if (i - start >= min_cnt) {
				uint32_t lq = 0, lr = 0, eq = 0, er = 0, sq = 0, sr = 0;
				mm_reg1_t *r;
				kv_pushp(mm_reg1_t, b->reg, &r);
				r->rid = rid, r->rev = rev, r->cnt = i - start, r->rep = 0;
				r->qs = ((uint32_t)b->mini.a[b->a[b->b[start]]>>32].y>>1) - (k - 1);
				r->qe = ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1) + 1;
				r->rs = rev? (uint32_t)b->a[b->b[i-1]] : (uint32_t)b->a[b->b[start]];
				r->re = rev? (uint32_t)b->a[b->b[start]] : (uint32_t)b->a[b->b[i-1]];
				r->rs -= k - 1;
				r->re += 1;
				r->mini_cnt = (b->a[b->b[i-1]]>>32) - (b->a[b->b[start]]>>32) + 1;
				for (j = start; j < i; ++j) { // count the number of times each minimizer is used
					int jj = b->a[b->b[j]]>>32;
					b->mini.a[jj].y += 1ULL<<32;
					kv_push(uint32_t, b->reg2mini, jj); // keep minimizer<=>reg mapping for derep
				}
				for (j = start; j < i; ++j) { // compute ->len
					uint32_t q = ((uint32_t)b->mini.a[b->a[b->b[j]]>>32].y>>1) - (k - 1);
					uint32_t r = (uint32_t)b->a[b->b[j]];
					r = !rev? r - (k - 1) : (0x80000000U - r);
					if (r > er) lr += er - sr, sr = r, er = sr + k;
					else er = r + k;
					if (q > eq) lq += eq - sq, sq = q, eq = sq + k;
					else eq = q + k;
				}
				lr += er - sr, lq += eq - sq;
				r->len = lr < lq? lr : lq;
			}
			start = i;
		}
	}
}

static inline void push_intv(mm128_v *intv, int start, int end) // used by get_reg() only
{
	mm128_t *p;
	if (intv->n > 0) {
		int last_start, last_end, min;
		p = &intv->a[intv->n-1];
		last_start = p->y, last_end = p->x + last_start;
		min = end - start < last_end - last_start? end - start : last_end - last_start;
		if (last_end > start && last_end - start > (min>>1) + (min>>2)) {
			p->x = end - last_start;
			return;
		}
	}
	kv_pushp(mm128_t, *intv, &p);
	p->x = end - start, p->y = start;
}

static void get_reg(mm_tbuf_t *b, int radius, int k, int min_cnt, int max_gap, int skip_derep)
{
	mm128_v *c = &b->coef;
	int i, start = 0;
	if (c->n < min_cnt) return;
	b->intv.n = 0;
	for (i = 1; i < c->n; ++i) { // identify all (possibly overlapping) clusters within _radius_
		if (c->a[i].x - c->a[start].x > radius) {
			if (i - start >= min_cnt) push_intv(&b->intv, start, i);
			for (++start; start < i && c->a[i].x - c->a[start].x > radius; ++start);
		}
	}
	if (i - start >= min_cnt) push_intv(&b->intv, start, i);
	radix_sort_128x(b->intv.a, b->intv.a + b->intv.n); // sort by the size of the cluster
	b->reg2mini.n = 0;
	for (i = b->intv.n - 1; i >= 0; --i) proc_intv(b, i, k, min_cnt, max_gap);
	if (!skip_derep) drop_rep(b, min_cnt);
}

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name)
{
	int j;

	b->mini.n = b->coef.n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, &b->mini);
	for (j = 0; j < b->mini.n; ++j) {
		int k, n;
		const uint64_t *r;
		int32_t qpos = (uint32_t)b->mini.a[j].y>>1, strand = b->mini.a[j].y&1;
		b->mini.a[j].y = b->mini.a[j].y<<32>>32; // clear the rid field
		r = mm_idx_get(mi, b->mini.a[j].x, &n);
		if (n > mi->max_occ) continue;
		for (k = 0; k < n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if ((opt->flag&MM_F_NO_SELF) && mi->name && strcmp(name, mi->name[r[k]>>32]) == 0 && rpos == qpos)
				continue;
			kv_pushp(mm128_t, b->coef, &p);
			if ((r[k]&1) == strand) { // forward strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
				p->y = (uint64_t)j << 32 | rpos;
			} else { // reverse strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos) | 1ULL<<63;
				p->y = (uint64_t)j << 32 | rpos;
			}
		}
	}
	radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
	b->reg.n = 0;
	get_reg(b, opt->radius, mi->k, opt->min_cnt, opt->max_gap, opt->flag&MM_F_WITH_REP);
	*n_regs = b->reg.n;
	return b->reg.a;
}

/**************************
 * Multi-threaded mapping *
 **************************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int batch_size, n_processed, n_threads;
	const mm_mapopt_t *opt;
	bseq_file_t *fp;
	const mm_idx_t *mi;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
	int *n_reg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	const mm_reg1_t *regs;
	int n_regs;

	regs = mm_map(step->p->mi, step->seq[i].l_seq, step->seq[i].seq, &n_regs, step->buf[tid], step->p->opt, step->seq[i].name);
	step->n_reg[i] = n_regs;
	if (n_regs > 0) {
		step->reg[i] = (mm_reg1_t*)malloc(n_regs * sizeof(mm_reg1_t));
		memcpy(step->reg[i], regs, n_regs * sizeof(mm_reg1_t));
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(s->n_seq, sizeof(int));
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t**));
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			for (j = 0; j < s->n_reg[i]; ++j) {
				mm_reg1_t *r = &s->reg[i][j];
				printf("%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
				if (mi->name) fputs(mi->name[r->rid], stdout);
				else printf("%d", r->rid + 1);
				printf("\t%d\t%d\t%d\t%d\t%d\t%.4f\n", mi->len[r->rid], r->rs, r->re, r->len, r->cnt, (double)r->cnt / r->mini_cnt);
			}
			free(s->reg[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
		}
		free(s->reg); free(s->n_reg); free(s->seq);
		free(s);
	}
    return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, int tbatch_size)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads, pl.batch_size = tbatch_size;
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	bseq_close(pl.fp);
	return 0;
}
