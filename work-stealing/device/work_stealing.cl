/**********************************************************************
 * Nadesh Ramanathan, Imperial College London
 *
 * File: work_stealing.cl 
 *
 * Additional Comments: distributed under a BSD license, see LICENSE.txt
 *
 **********************************************************************/
typedef struct {
    //Optimisation to reduce Altera BRAM replication
    //    uint tail;
    //    uint head;
    //    bool finish;
    task_type stack[STACK_SIZE];
} deque_type;

void cl_push(task_type stack, 
        local uint * restrict tail, 
        local task_type (* restrict deque)[STACK_SIZE/4]
        ){

    uint tmp = *tail;
    short unsigned int col = tmp >> SHIFT; 
    deque[col][tmp]=stack;
    tmp++;
    *tail = tmp;
}

bool cl_pop(task_type * task_record, 
        local uint * restrict tail,
        local uint * restrict head,
        local task_type (* restrict deque)[STACK_SIZE/4]
        ){

#ifndef SYN
    uint par = get_local_id(0);
#endif

    uint oldTail = *tail; 
    uint oldHead = *head;
    short unsigned int col = oldTail >> SHIFT; 

    if(oldTail==0){
#ifdef VERBOSE
#ifndef SYN
        printf("%d: Stack Empty\n",par);
#endif
#endif
        return false;
    }

    uint oldHead_index = oldHead & 0x0000FFFF;
    uint oldHead_ctr = oldHead >> 16;
    uint newHead = (oldHead_ctr+1) << 16;
    uint tmp = oldTail - 1;

    *task_record = deque[col][tmp]; 

    *tail = (tmp > oldHead_index) ? tmp : 0;

    if(tmp > oldHead_index){
        return true;
    }

    if(tmp == oldHead_index){
#ifdef VERBOSE
#ifndef SYN
        printf("%d: Non Empty One element T: %d, H: %d, Ctr: %d, Newhead %x\n",par,tmp,oldHead_index, oldHead_ctr, newHead);
#endif
#endif

        if(oldHead == atomic_cmpxchg(head,oldHead,newHead)){
#ifdef VERBOSE
#ifndef SYN
            printf("%d: One-element pop successful\n",par);
#endif
#endif
            return true;
        }

    }

#ifdef VERBOSE
#ifndef SYN
    printf("%d: Stack Empty 1\n",par);
#endif
#endif
    *head = newHead;
    return false;

}

bool cl_steal(task_type * task_record, 
        local uint * restrict tail,
        local uint * restrict head,
        local task_type (* restrict deque)[STACK_SIZE/4]
        ){

#ifndef SYN
    uint par = get_local_id(0);
#endif
    uint oldHead = *head;
    uint oldTail = *tail;
    short unsigned int col = oldTail >> SHIFT; 
    uint oldHead_index = oldHead & 0x0000FFFF;
    uint oldHead_ctr = oldHead >> 16;
    uint newHead = oldHead;
    newHead++;

    *task_record = deque[col][oldHead_index];

    if( oldTail <= oldHead_index){
#ifdef VERBOSE
#ifndef SYN
        printf("Steal failure. Stack %d Empty\n",par);
#endif
#endif
        return false;
    }

    if(oldHead == atomic_cmpxchg(head,oldHead,newHead)){
#ifdef VERBOSE
#ifndef SYN
        printf("Stealing successful\n");
#endif
#endif
        return true;
    }

#ifdef VERBOSE
#ifndef SYN
    printf("Stealing unsuccessful although queue not empty, H: %d\n",*head);
#endif
#endif
    return false;

}

bool cl_get(uint local_id, 
        uint steal_id, 
        uint * next_steal_id, 
        uint * return_id,
        bool * try_to_steal, 
        task_type * task_record, 
        local uint * restrict tail,
        local uint * restrict head,
        local task_type (* restrict deque)[4][STACK_SIZE/4]
        ){


#ifdef CROSSBAR
#if P > 2
#ifdef VERBOSE
#ifndef SYN 
    printf("%d: Initial p %d ", local_id, *next_steal_id);
#endif
#endif
    *next_steal_id = (*next_steal_id+1)%P;
    if(*next_steal_id==local_id){ 
        *next_steal_id=(*next_steal_id+1)%P;
    }
#ifdef VERBOSE
#ifndef SYN 
    printf("Final p %d\n",*next_steal_id);
#endif
#endif
#endif
#endif

    *return_id = local_id;
    uint success = false;
    success = cl_pop(task_record, &tail[local_id], &head[local_id], deque[local_id]);

    // return success; 

    if(success) {
        return true;
    }
    *try_to_steal = true;

#ifndef SYN
    printf("Stack %d empty -> Trying to steal from %d Tail of stack %d: T: %d H: %d\n", local_id, steal_id, steal_id, tail[steal_id], (head[steal_id] &0x0000ffff)  );
#endif

    success = cl_steal(task_record, &tail[steal_id], &head[steal_id], deque[steal_id]);

    if(success) *return_id = steal_id;

    return success;
}
