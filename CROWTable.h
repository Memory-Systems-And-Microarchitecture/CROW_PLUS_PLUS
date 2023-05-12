#ifndef __H_CROW_TABLE_H__
#define __H_CROW_TABLE_H__

#include <unordered_map>
#include <list>
#include <cstdlib>
#include <ctime>
#include "Statistics.h"
#include <fstream>
#include <iostream>


using namespace std;

namespace ramulator {

    class CROWEntry {
        public:
            ulong row_addr;
            bool valid;
            bool FR; //Force Restoration
            int hit_count;
            long total_hits; // used only for statistic collection
            bool is_to_remap_weak_row;
            bool    make_it_static;
            CROWEntry() {
                row_addr = 0;
                valid = false;
                FR = false;
                hit_count = 0;
                total_hits = 0;
                is_to_remap_weak_row = false;
                make_it_static = false;
            }
    };
  
    
    // crow++ starts
    class CROWLifetimeEntry {
    public:
        ulong row_addr;
        //uint hit_count;
        int entry_into_cache_count;
        int bank;
        int hit_count;
        CROWLifetimeEntry(){
            row_addr=0;
            bank=0;
            hit_count = 0;
            //hit_count = 0;
            entry_into_cache_count=0;
        }
    };
    // crow++ ends
 
    template <typename T>
    class CROWTable {


    protected:
        ScalarStat crow_evict_with_zero_hits;
        ScalarStat crow_evict_with_one_hit;
        ScalarStat crow_evict_with_two_hits;
        ScalarStat crow_evict_with_5_or_less_hits;
        ScalarStat crow_evict_with_6_or_more_hits;
         
    public:

        bool is_DDR4 = false, is_LPDDR4 = false;
        bool enable_crow_plus_plus;
        // crow++
             //int track_total_entries = 0;
             //track_total_entries = (spec->org_entry.count[int(T::Level::Bank)]) * num_SAs * num_track_rows;
             //cout<<"total_entries inside new table:"<<track_total_entries<<endl;
             CROWLifetimeEntry* crow_track_evict_entry = new CROWLifetimeEntry [8*128*8];
        int temp_track_add_entry = 0;
       
        // crow++
        CROWTable(const T* spec, const int crow_id, const uint num_SAs, 
                    const uint num_copy_rows, const uint num_weak_rows, 
                    const uint crow_evict_hit_thresh,
                    const uint crow_half_life, const float to_mru_frac,
                    const uint num_grouped_SAs) : 
                        spec(spec), crow_id(crow_id), num_SAs(num_SAs), num_copy_rows(num_copy_rows),
                        num_weak_rows(num_weak_rows), to_mru_frac(to_mru_frac), num_grouped_SAs(num_grouped_SAs){
                    // num_grouped_SAs indicates how many consecutively
                    // addressed subarrays (SAs) will share CROW table entries. This is
                    // implemented to optimize the storage requirement of
                    // the CROW table. num_grouped_SAs is 1 by default,
                    // which means each SA will have its own num_copy_rows number
                    // of entries.

            assert((num_weak_rows <= num_copy_rows) && "Cannot have more weak rows than copy rows! Need to fallback to the default refresh rate in case there are too many weak rows.");

             
             

            int num_levels = int(T::Level::Bank) - int(T::Level::Rank) + 1;

            sizes = new int[num_levels + 1]; // +1 for subarrays


            num_entries = num_SAs * num_copy_rows;
            for(int i = 0; i < num_levels; i++) {
                sizes[i] = spec->org_entry.count[1 + i]; // +1 to skip channel
                num_entries *= spec->org_entry.count[1 + i]; 
            }

            num_entries /= num_grouped_SAs;

            // crow++
            //  for(int i = 0; i < (spec->org_entry.count[int(T::Level::Bank)]); ++i) {
            // crow_track_evict_entry[i] = new CROWLifetimeEntry*[num_SAs];
            // for (int j = 0; j < num_SAs ; ++j) {
            // crow_track_evict_entry[i][j] = new CROWLifetimeEntry[num_track_rows];
            //     }
            // }
            // crow++

            entries = new CROWEntry[num_entries];
            lru_lists = new list<CROWEntry*>[num_entries/num_copy_rows];
            pointer_maps = new unordered_map<CROWEntry*, 
                         list<CROWEntry*>::iterator>[num_entries/num_copy_rows];

            cur_accesses = new int[num_entries/num_copy_rows];
            fill(cur_accesses, cur_accesses + (num_entries/num_copy_rows), 0);
            hit_count_half_life = crow_half_life;

            srand(static_cast<unsigned>(0));

            next_copy_row_id = new uint[num_entries/num_copy_rows];
            for(int i = 0; i < (num_entries/num_copy_rows); i++) {
                next_copy_row_id[i] = 0;
            }

            is_DDR4 = spec->standard_name == "DDR4";
            is_LPDDR4 = spec->standard_name == "LPDDR4";

            if(is_DDR4 || is_LPDDR4)
                SA_size = spec->org_entry.count[int(T::Level::Row)]/num_SAs;
            else // for SALP
                SA_size = spec->org_entry.count[int(T::Level::Row) - 1];

            sizes[num_levels] = num_copy_rows;
            //cout<<"sizes_num_levels: "<<sizes[num_levels]<<endl;
            sizes[num_levels-1] = sizes[num_levels] * (num_SAs/num_grouped_SAs);
            //cout<<"sizes_num_levels-1: "<<sizes[num_levels-1]<<endl;
            for(int i = (num_levels - 2); i >= 0; i--) {
                sizes[i] = spec->org_entry.count[i + 2] * sizes[i + 1];
                //cout<<"sizes["<<i<<"]: "<<sizes[i]<<endl;
            }

            do_weak_row_mapping();

            //for(int i = 0; i <= num_levels; i++)
                //printf("sizes[%d] = %d \n", i, sizes[i]);

            crow_evict_with_zero_hits
                .name("crow_evict_with_zero_hits_"+to_string(crow_id))
                .desc("Number of copy rows evicted without being hit")
                .precision(0)
                ;

            crow_evict_with_one_hit
                .name("crow_evict_with_one_hit_"+to_string(crow_id))
                .desc("Number of copy rows evicted with a single hit")
                .precision(0)
                ;
            crow_evict_with_two_hits
                .name("crow_evict_with_two_hits_"+to_string(crow_id))
                .desc("Number of copy rows evicted with two hits")
                .precision(0)
                ;
            crow_evict_with_5_or_less_hits
                .name("crow_evict_with_5_or_less_hits_"+to_string(crow_id))
                .desc("Number of copy rows evicted with 5 hits or less")
                .precision(0)
                ;
            crow_evict_with_6_or_more_hits
                .name("crow_evict_with_6_or_more_hits_"+to_string(crow_id))
                .desc("Number of copy rows evicted with 6 or more hits")
                .precision(0)
                ;


        }

        virtual ~CROWTable() {
            delete[] entries;
            delete[] next_copy_row_id;
            delete[] sizes;

            delete[] lru_lists;
            delete[] pointer_maps;
            delete[] cur_accesses;
        }

        bool access(const vector<int>& addr_vec, const bool move_to_LRU = false) {

            int ind = (calc_entries_offset(addr_vec)/num_copy_rows);
            int& cur_access = cur_accesses[ind];
            cur_access++;

            if(cur_access == hit_count_half_life) {
                cur_access = 0;

                // traverse all corresponding entries and halve their hit counts
                auto& cur_lru_list = lru_lists[ind];

                for(auto it = cur_lru_list.begin(); it != cur_lru_list.end(); it++) {
                    (*it)->hit_count >>= 1;
                } 
            }

            if(is_hit(addr_vec)) {
                // promote to most recently accessed
                CROWEntry* cur_entry = get_hit_entry(addr_vec);

                cur_entry->total_hits++;
                if(cur_entry->hit_count < MAX_HIT_COUNT)
                    cur_entry->hit_count++;

                auto& cur_lru_list = lru_lists[ind];
                auto& cur_pointer_map = pointer_maps[ind];

                cur_lru_list.erase(cur_pointer_map[cur_entry]);
                if(!move_to_LRU) {
                    cur_lru_list.push_front(cur_entry);
                    cur_pointer_map[cur_entry] = cur_lru_list.begin();
                } else {
                    cur_lru_list.push_back(cur_entry);
                    cur_pointer_map[cur_entry] = prev(cur_lru_list.end());
                }

                return true;
            }
            else {
                assert(!move_to_LRU && "Error: Move to LRU should not be set on a miss!");
                // do not insert to lru_list on miss
                // insert only when add_entry() is called
                
                // decrease the hit counter of the LRU entry
                //CROWEntry* LRU_entry = get_LRU_entry(addr_vec);
                //if(LRU_entry->hit_count > 0)
                //    LRU_entry->hit_count--;
            }

            return false;

        }

        void make_LRU(const vector<int>& addr_vec, CROWEntry* crow_entry) {
            int ind = (calc_entries_offset(addr_vec)/num_copy_rows);
            auto& cur_lru_list = lru_lists[ind];
            auto& cur_pointer_map = pointer_maps[ind];

            cur_lru_list.erase(cur_pointer_map[crow_entry]);
            cur_lru_list.push_back(crow_entry);
            cur_pointer_map[crow_entry] = prev(cur_lru_list.end());
        }

        // this function does not update the LRU state
        bool is_hit(const vector<int>& addr_vec) {
            
            if ((addr_vec[int(T::Level::Bank)] == 6) and (addr_vec[int(T::Level::Row)] / SA_size == 94)){
                //cout<<"found Bank 6 and SA = 94"<<endl;
                //cout<<"Problem here please check"<<endl;
                //cout<<"we need: "<<addr_vec[int(T::Level::Row)]<<endl;
                bool track_enter = false;
                if(get_hit_entry(addr_vec) == nullptr){
                    for(uint i = 0; i < num_copy_rows; i++) {
                      // cout<<"row_addr: "<<get_entry(addr_vec, i)->row_addr<<endl;
                    }
                    track_enter = true;
                }
                if(track_enter){
                    //cout<<"not found"<<endl;
                } else{
                    //cout<<"found"<<endl;
                }
            }


            if(get_hit_entry(addr_vec) != nullptr){
                if(num_weak_rows == num_copy_rows)
                    assert(false && "WARNING: We should not get hit when all copy rows are allocated for weak rows.");
                return true;
            }

            return false;
        }

        CROWEntry* add_entry(const vector<int>& addr_vec, const bool FR) {
            temp_track_add_entry+=1;
           // cout<<"adding a new row"<<endl;
            int lru_ind = calc_entries_offset(addr_vec)/num_copy_rows;
            auto& cur_lru_list = lru_lists[lru_ind];
            auto& cur_pointer_map = pointer_maps[lru_ind];
            CROWEntry* cur_entry = nullptr;

            int freeLoc = free_loc(addr_vec);
            if (freeLoc == -1) {
                cur_entry = get_LRU_entry(addr_vec);
                
                //crow++
                   int index_entries = 0;
                   for(;index_entries<(8*128*8);index_entries++){
                      if(crow_track_evict_entry[index_entries].bank == addr_vec[int(T::Level::Bank)] and crow_track_evict_entry[index_entries].row_addr == cur_entry->row_addr){
                              // cout<<"adding count"<<endl;
                              // cout<<"row_addr: "<<cur_entry->row_addr<<endl;
                              // cout<<"bank: "<<addr_vec[int(T::Level::Bank)]<<endl;
                               crow_track_evict_entry[index_entries].hit_count = cur_entry->total_hits;
                              // cout<<crow_track_evict_entry[index_entries].hit_count<<endl;
                               if (crow_track_evict_entry[index_entries].hit_count>5 and crow_track_evict_entry[index_entries].entry_into_cache_count>2){
                                      cur_entry->make_it_static = true;
                                    
                               }

                               break;
                      }
                   }

                //crow++

                assert(!cur_entry->FR && "The discarded entry should not require full restoration!");
                
                if(cur_entry->total_hits == 0){
                    crow_evict_with_zero_hits++;
                    crow_evict_with_5_or_less_hits++;
                }
                else if(cur_entry->total_hits == 1){
                    crow_evict_with_one_hit++;
                    crow_evict_with_5_or_less_hits++;
                }
                else if(cur_entry->total_hits == 2){
                    crow_evict_with_two_hits++;
                    crow_evict_with_5_or_less_hits++;
                }
                else if(cur_entry->total_hits <= 5)
                    crow_evict_with_5_or_less_hits++;
                else
                    crow_evict_with_6_or_more_hits++;
                //cout<<"Evicting entry"<<endl;
                store_each_evict_entry(cur_entry, addr_vec[int(T::Level::Bank)]);
                cur_entry->row_addr = addr_vec[int(T::Level::Row)];
                cur_entry->valid = true;
                cur_entry->FR = FR;
                cur_entry->hit_count = 0;
                cur_entry->total_hits = 0;

                update_next_copy_row_id(addr_vec);
                store_each_add_entry(cur_entry,addr_vec[int(T::Level::Bank)]);
                // remove the LRU entry from the list
                assert(cur_lru_list.size() <= num_copy_rows);

                                
                cur_lru_list.pop_back();
                cur_pointer_map.erase(cur_entry);
            } 
            else {
                cur_entry = get_entry(addr_vec, uint(freeLoc));
                cur_entry->row_addr = addr_vec[int(T::Level::Row)];
                cur_entry->valid = true;
                cur_entry->FR = FR;
                cur_entry->hit_count = 0;
                cur_entry->total_hits = 0;
                
                // As there is a free location, cur_lru_list should not be
                // full
                store_each_add_entry(cur_entry,addr_vec[int(T::Level::Bank)]);
                assert(cur_lru_list.size() < num_copy_rows && 
                        "There is a free copy row, the LRU table should not be full!"); 
            }

            // crow++

             //check if already added 
             int index_entries = 0;
             int threshold_least_count = 100;
             int least_count_index;
             int temp_size = 8*128*8;
             for(;index_entries<temp_size;index_entries++){
                if(crow_track_evict_entry[index_entries].bank == addr_vec[int(T::Level::Bank)] and crow_track_evict_entry[index_entries].row_addr == cur_entry->row_addr){
                    //already present row coming again
                    crow_track_evict_entry[index_entries].entry_into_cache_count += 1;
                    switch(crow_track_evict_entry[index_entries].entry_into_cache_count){
                        case 1:
                           if (enable_crow_plus_plus)
                              to_mru_frac = 0.5;
                           break;
                        case 2:
                            if (enable_crow_plus_plus)
                               to_mru_frac = 0.5;
                            break;
                        case 3:
                            if (enable_crow_plus_plus)
                               to_mru_frac = 0.7;
                        case 4:
                            if(enable_crow_plus_plus)
                               to_mru_frac = 0.8;
                        case 5:
                            //cout<<"Goal Hit Make it Static"<<endl;
                            to_mru_frac = 0.9;
                            //cout<<"row_addr: "<<cur_entry->row_addr<<" hits: "<<crow_track_evict_entry[index_entries].entry_into_cache_count<<" hit_count: "<<crow_track_evict_entry[index_entries].hit_count<<endl;
                    }
                    
                    
                    break;
                }
                //check least threshold count for replacement inside this table
                    if (crow_track_evict_entry[index_entries].entry_into_cache_count < threshold_least_count){
                        threshold_least_count = crow_track_evict_entry[index_entries].entry_into_cache_count;
                        least_count_index = index_entries;
                    } 
             }
            // check for if not address found in previous case means a fresh entry
            if (index_entries>=temp_size){
                 //cout<<"Address Not Found In Cache: "<<endl;
                 //cout<<"row_addr: "<<crow_track_evict_entry[least_count_index].row_addr<<" count: "<<crow_track_evict_entry[least_count_index].hit_count<<" entry_into_cache_count: "<<crow_track_evict_entry[least_count_index].entry_into_cache_count<<endl;
                 crow_track_evict_entry[least_count_index].row_addr = cur_entry->row_addr;
                 crow_track_evict_entry[least_count_index].bank =  addr_vec[int(T::Level::Bank)];
                 crow_track_evict_entry[least_count_index].entry_into_cache_count = 1;
            } else{
                //cout<<"Address already Present in cache: "<<endl;
                //cout<<"row_addr: "<<crow_track_evict_entry[index_entries].row_addr<<" count: "<<crow_track_evict_entry[index_entries].hit_count<<" entry_into_cache_count: "<<crow_track_evict_entry[index_entries].entry_into_cache_count<<endl;
                if (crow_track_evict_entry[index_entries].hit_count >= 2 and crow_track_evict_entry[index_entries].entry_into_cache_count >= 2 and enable_crow_plus_plus == true){
                             //cout<<"hit"<<endl;
                             to_mru_frac = 1;
                }
                // crow_track_evict_entry[least_count_index].row_addr = cur_entry->row_addr;
                // crow_track_evict_entry[least_count_index].bank =  addr_vec[int(T::Level::Bank)];
                // crow_track_evict_entry[least_count_index].entry_into_cache_count = 1;
            }

            //  int ii=0;
            //  for(;ii<num_track_rows;ii++){
            //     //cout<<"rank: "<<addr_vec[int(T::Level::Rank)]<<endl;
            //     //cout<<"bank: "<<addr_vec[int(T::Level::Bank)]<<endl;
            //     //cout<<"row: "<<addr_vec[int(T::Level::Row)]<<endl;
            //     //cout<<"subarray: "<<addr_vec[int(T::Level::Row)] / SA_size<<endl;
            //     //cout<<"matched address: "<<crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].row_addr<<endl;
            //     //cout<<"required address: "<<cur_entry->row_addr<<endl;
               
            //     if(crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].row_addr == cur_entry->row_addr){
            //         crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count += 1;
            //        // cout<<"address inside our table: "<<crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].row_addr<<endl;
            //        // cout<<"adress required: "<<cur_entry->row_addr<<endl;
            //         //cout<<"Added one"<<crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count<<endl;
            //         switch (crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count) {
            //         case 1:
            //             //cout<<"coming second time"<<endl;
            //             if (enable_crow_plus_plus)
            //                 to_mru_frac = 0.6; // 0.1
            //             break;
            //         case 2:
            //             //cout<<"Coming third time"<<endl;
            //             if (enable_crow_plus_plus)
            //                 to_mru_frac = 0.7;  // 0.5
            //             break;
            //         case 3:
            //             //cout<<"Coming fourth time"<<endl;
            //             if (enable_crow_plus_plus)
            //                 to_mru_frac = 0.8;  // 0.7
            //             break;
            //         case 4:
            //             //cout<<"Coming fifth time"<<endl;
            //             if (enable_crow_plus_plus)
            //                 to_mru_frac = 0.9;  // 0.8
            //             break;
            //         case 5:
            //             //cout<<"Coming sixth time"<<endl;
            //             if (enable_crow_plus_plus)
            //                 to_mru_frac = 1.0;
            //             break;
            //         default:
            //             to_mru_frac = 0.0;
            //             break;
            //         }
            //         if (crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count == 10){
            //             cout<<"Goal Hit"<<endl;
            //             //cur_entry->reserved_lifetime = true;
            //         }
            //         break;
            //     }
            //  }
            //  if (ii>=num_track_rows){
            //     int ii=0;
            //     int least_val = 100;
            //     int index_store = 0;
            //     for(;ii<=num_track_rows;ii++){
            //         if (crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count < least_val){
            //             least_val = crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][ii].entry_into_cache_count;
            //             index_store = ii;
            //         }
            //     }
            //     //cout<<"Before Adding a row: "<<crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][index_store].row_addr<<endl;
            //     crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][index_store].row_addr = cur_entry->row_addr;
            //    // cout<<"After Adding a row: "<<crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][index_store].row_addr<<endl;
            //     crow_track_evict_entry[addr_vec[int(T::Level::Bank)]][(addr_vec[int(T::Level::Row)] / SA_size)][index_store].entry_into_cache_count = 1;
            //     to_mru_frac = 0.0;
            //     store_entry_to_file(addr_vec[int(T::Level::Bank)],(addr_vec[int(T::Level::Row)] / SA_size), cur_entry->row_addr);
            //     for(int t1 = 0; t1 < (spec->org_entry.count[int(T::Level::Bank)]); ++t1) {
            //         //cout<<"Printing Details of Bank: "<<t1<<endl;
            //     for (int t2 = 0; t2 < num_SAs ; ++t2) {
            //         //cout<<"Printing Details of Subarray: "<<t2<<endl;
            //         for(int t3=0;t3<num_track_rows;t3++){
            //           //  cout<<"row_address: "<<crow_track_evict_entry[t1][t2][t3].row_addr<<" Count: "<<crow_track_evict_entry[t1][t2][t3].entry_into_cache_count<<endl;
            //         }
            //     }
            //     //cout<<"Bank -->"<<" "<<i<<" created"<<endl;
            //     //cout<<"Number of SA -->"<<" "<<num_SAs<<endl;
            //     //cout<<"Number of Rows Per SA -->"<<" "<<num_track_rows<<endl;
            //     }
            //  }

            // crow++

            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            if(r >= to_mru_frac){
                cur_lru_list.push_back(cur_entry);
                cur_pointer_map[cur_entry] = prev(cur_lru_list.end());
            } else {
                //cout<<"storing at MRU"<<endl;
                store_at_mru(addr_vec[int(T::Level::Bank)], cur_entry->row_addr, (addr_vec[int(T::Level::Row)] / SA_size), crow_track_evict_entry[index_entries].hit_count);
                //cout<<"Inserting at MRU, row_addr: "<<cur_entry->row_addr<<endl;
                cur_lru_list.push_front(cur_entry);
                cur_pointer_map[cur_entry] = cur_lru_list.begin();
            }

            return cur_entry;
        }

        CROWEntry* get_entry(const vector<int>& addr_vec, const uint copy_row_id) {
            int offset = calc_entries_offset(addr_vec);    
            return &entries[offset + copy_row_id];
        }

        CROWEntry* get_hit_entry(const vector<int>& addr_vec){
            int offset = calc_entries_offset(addr_vec);

            int row_addr = addr_vec[int(T::Level::Row)];
            for (uint i = 0; i < num_copy_rows; i++) {
                if (entries[offset + i].row_addr == row_addr && 
                        entries[offset + i].valid && !entries[offset + i].is_to_remap_weak_row){
                    store_each_hit_entry(&entries[offset+i], addr_vec[int(T::Level::Bank)], (addr_vec[int(T::Level::Row)] / SA_size));
                    return &entries[offset + i];
                }
            }

            return nullptr;
        }

        bool set_FR(const vector<int>& addr_vec, const bool FR) {
            CROWEntry* cur_entry = get_hit_entry(addr_vec);
        
            if(cur_entry == nullptr){
                assert(false && "We should always get a hit.");
                return false;
            }
        
            cur_entry->FR = FR;
        
            return true;
        }

        CROWEntry* get_LRU_entry(const vector<int>& addr_vec, int crow_evict_threshold = 0){
            int lru_ind = calc_entries_offset(addr_vec)/num_copy_rows;
            auto& cur_lru_list = lru_lists[lru_ind];
            
            if (enable_crow_plus_plus){
            auto itr = cur_lru_list.end();
            auto itr_return = cur_lru_list.back();
            int count = 100;
            do{
              itr--;
              if ((*itr)->hit_count < count){
                count = (*itr)->hit_count;
                itr_return = (*itr);
              }
              //std::cout<<"row_addr: "<<(*itr)->row_addr<<" hit_count: "<<(*itr)->hit_count<<endl;
            } while(itr != cur_lru_list.begin());
            //cout<<"LRU item: "<<itr_return->hit_count<<endl;
            store_LRU(itr_return->row_addr,itr_return->hit_count);
            //updated changes
           
            //updated changes
            return itr_return;
           }
           if(crow_evict_threshold == 0){
                //cout<<"LRU item: "<<cur_lru_list.back()->hit_count<<endl;
                store_LRU(cur_lru_list.back()->row_addr,cur_lru_list.back()->hit_count);
                return cur_lru_list.back();
           }

            auto it = cur_lru_list.end();
            do {
                it--;
                if((*it)->hit_count <= crow_evict_threshold)
                   return *it;
            } while(it != cur_lru_list.begin());
            
            return nullptr;
        }

        CROWEntry* get_discarding_entry(const vector<int>& addr_vec){
            return get_entry(addr_vec, get_next_copy_row_id(addr_vec));
        }

        bool get_discarding_FR(const vector<int>& addr_vec) {
            return get_entry(addr_vec, get_next_copy_row_id(addr_vec))->FR;
        }

        ulong get_discarding_row_addr(const vector<int>& addr_vec) {
            return get_entry(addr_vec, get_next_copy_row_id(addr_vec))->row_addr;
        }

        int get_discarding_copy_row_id(const vector<int>& addr_vec) {
            return get_next_copy_row_id(addr_vec);
        }

        bool is_full(const vector<int>& addr_vec) {
           return (free_loc(addr_vec) == -1); 
        }

        void invalidate(const vector<int>& addr_vec) {
            CROWEntry* entry = get_hit_entry(addr_vec);
            assert(!entry->is_to_remap_weak_row && "A remapped weak row should not be invalidated!");
            int lru_ind = calc_entries_offset(addr_vec)/num_copy_rows;
            auto& cur_lru_list = lru_lists[lru_ind];
            auto& cur_pointer_map = pointer_maps[lru_ind];

            cur_lru_list.erase(cur_pointer_map[entry]);
            cur_pointer_map.erase(entry);

            entry->valid = false;

        }
        
        int find_not_FR(const vector<int>& addr_vec) {
            int offset = calc_entries_offset(addr_vec);

            for(int i = 0; i < num_copy_rows; i++) {
                if(!entries[offset + i].FR)
                    return i;
            }

            return -1;
        }
        
        void print() {
        
            for(int i = 0; i < num_entries; i++) {
                CROWEntry& ce = entries[i];
                printf("CROWTable: %lu %d %d \n", ce.row_addr, ce.valid, ce.FR);
            }
        }

        // crow++
       
        
        void initialize_trace_file(){
            outfile_trace_add_entry.open("crowtable_trace_add_entry.txt");
            outfile_trace_evict_entry.open("crowtable_trace_evict_entry.txt");
            outfile_trace_hit_entry.open("crowtable_trace_hit_entry.txt");
            outfile_trace_makeLRU_entry.open("crowtable_trace_makeLRU_entry.txt");
            outfile_trace_makeEntryLifetime.open("crowtable_trace_makeEntryLifetime.txt");
            outfile_store_proposed_table.open("trace_of_proposed_table.txt");
            outfile_store_insert_at_mru.open("trace_insert_at_mru.txt");
            if (enable_crow_plus_plus)
               outfile_store_LRU.open("trace_LRU_crow_plus_plus.txt");
            else
               outfile_store_LRU.open("trace_LRU_crow.txt");
            first_entry_add = false;
        }

        void finish_all(){
            outfile.close();
            outfile_trace_add_entry.close();
            outfile_trace_evict_entry.close();
            outfile_trace_hit_entry.close();
            outfile_store_proposed_table.close();
            outfile_store_insert_at_mru.close();
            outfile_store_LRU.close();
        }
        
         void store_entry_to_file(int bank, int sa, int row_addr){
            if(first_entry_add){
              initialize_trace_file();        
            }
            outfile_store_proposed_table<<bank<<" "<<sa<<" "<<row_addr<<endl;
        }
        
        void store_LRU(long row_addr, int hit_count){
           if(first_entry_add){
            initialize_trace_file();
           }
           outfile_store_LRU<<"row_addr: "<<row_addr<<" hit_count: "<<hit_count<<endl;
        }

        void store_into_file() {
            ofstream outfile;
            outfile.open("crowtable_trace.txt");
            outfile<<"Crow table"<<endl;
            for(int i = 0; i < num_entries; i++) {
                CROWEntry& ce = entries[i];
                //printf("CROWTable: %lu %d %d \n", ce.row_addr, ce.valid, ce.FR);
                outfile<<"row_addr: "<<ce.row_addr<<" valid: "<<ce.valid<<" force_restoration: "<<ce.FR<<endl;
            }
        }

        void store_each_add_entry(CROWEntry* ce, int bank){
            if (first_entry_add){
                initialize_trace_file();
            }
            outfile_trace_add_entry<<"BanK: "<<bank<<" SA_addr: "<<(ce->row_addr/SA_size)<<" row_addr: "<<ce->row_addr<<" hits: "<<ce->hit_count<<endl;
        }
         
        void store_each_evict_entry(CROWEntry* ce, int bank){
            if (first_entry_add){
                initialize_trace_file();
            }
            //cout<<"Inside store evict entry"<<endl;
            outfile_trace_evict_entry<<"Bank: "<<bank<<" SA_addr: "<<(ce->row_addr/SA_size)<<" row_addr: "<<ce->row_addr<<" hits: "<<ce->hit_count<<endl;
        }



        void store_each_hit_entry(CROWEntry* ce, int bank, int sa){
            if (first_entry_add){
                initialize_trace_file();
            }
            outfile_trace_hit_entry<<"Bank: "<<bank<<" SA_addr: "<<(ce->row_addr / SA_size)<<" row_addr: "<<ce->row_addr<<" hits: "<<ce->hit_count<<endl;
        }

        void store_makeLRU_entry(CROWEntry* ce){
            if(first_entry_add){
                initialize_trace_file();
            }
            outfile_trace_makeLRU_entry<<"row_addr: "<<ce->row_addr<<" valid: "<<ce->valid<<" hit_count: "<<ce->hit_count<<" total_hits: "<<ce->total_hits<<" force_restoration: "<<ce->FR<<endl;
        }
        void store_at_mru(int bank, int row_addr, int subarrays_addr, int hit_count){
            if(first_entry_add){
                initialize_trace_file();
            }
            outfile_store_insert_at_mru<<"Bank: "<<bank<<" SA_addr: "<<subarrays_addr<<" row_addr: "<<row_addr<<" hit_count: "<<hit_count<<endl; 
        }

        //crow++

    private:
        const T* spec;
        int crow_id;
        CROWEntry* entries;
        int num_entries = 0;
        int* sizes;
        uint* next_copy_row_id;
        uint num_SAs = 1;
        uint num_copy_rows;
        uint num_weak_rows;
        int SA_size = 0;
        
        // Crow++
        uint num_track_rows = 8;
        ofstream outfile;
        bool first_entry_add = true;
        ofstream outfile_trace_add_entry;
        ofstream outfile_trace_evict_entry;
        ofstream outfile_trace_hit_entry;
        ofstream outfile_trace_makeLRU_entry;
        ofstream outfile_trace_makeEntryLifetime;
        ofstream outfile_store_proposed_table;
        ofstream outfile_store_insert_at_mru;
        ofstream outfile_store_LRU;
       
        //

        float to_mru_frac = 0.6f; //by default, always insert a new entry to LRU position
        uint num_grouped_SAs = 1;

        int hit_count_half_life = 1000; // in number of CROWTable accesses
        const int MAX_HIT_COUNT = 32;
        int* cur_accesses = nullptr;

        // LRU replacement policy
        list<CROWEntry*>* lru_lists; //one for each subarray
        unordered_map<CROWEntry*, list<CROWEntry*>::iterator>* pointer_maps;


        void update_next_copy_row_id(const vector<int>& addr_vec) {
            uint cur = get_next_copy_row_id(addr_vec);

            int offset = calc_next_copy_row_id_offset(addr_vec);
            next_copy_row_id[offset] = (cur + 1) % num_copy_rows;
        }

        uint get_next_copy_row_id(const vector<int>& addr_vec) {
            int offset = calc_next_copy_row_id_offset(addr_vec);
            return next_copy_row_id[offset]; 
        }

        int free_loc(const vector<int>& addr_vec) {
            if ((addr_vec[int(T::Level::Bank)] == 6) and (addr_vec[int(T::Level::Row)] / SA_size == 94)){
               // cout<<"found Bank 6 and SA = 94"<<endl;
                for(uint i = 0; i < num_copy_rows; i++) {
                   // cout<<(get_entry(addr_vec, i)->row_addr / 128)<<" valid: "<<get_entry(addr_vec, i)->valid<<endl;
                }
            }
            for(uint i = 0; i < num_copy_rows; i++) {
                if(!(get_entry(addr_vec, i)->valid))
                    return i;
            }

            return -1;
        }

        int calc_entries_offset(const vector<int>& addr_vec) {
            int offset = 0;

            // skipping channel
            for(int l = 1; l <= int(T::Level::Bank); l++) {
                offset += sizes[l-1]*addr_vec[l];
            }

            int SA_id;
            if(is_DDR4 || is_LPDDR4)
                SA_id = addr_vec[int(T::Level::Row)]/SA_size;
            else // for SALP
                SA_id = addr_vec[int(T::Level::Row) - 1];
            offset += sizes[int(T::Level::Bank)]*(SA_id/num_grouped_SAs);
            //cout<<"Bank: "<<addr_vec[int(T::Level::Bank)]<<" Row: "<<addr_vec[int(T::Level::Row)]<<" SA_Id: "<<(addr_vec[int(T::Level::Row)] / SA_size)<<" offset: "<<offset<<endl;
            return offset;       
        }


        int calc_next_copy_row_id_offset(const vector<int>& addr_vec) {
           int offset = 0;
           // skipping channel
           for(int l = 1; l <= int(T::Level::Bank); l++) {
               offset += (sizes[l-1]*addr_vec[l])/num_copy_rows;
           }
       
           int SA_id;
           if(is_DDR4 || is_LPDDR4)
               SA_id = addr_vec[int(T::Level::Row)]/SA_size;
           else // for SALP
               SA_id = addr_vec[int(T::Level::Row) - 1];
           offset += (SA_id/num_grouped_SAs);
       
           return offset;
        }

        void do_weak_row_mapping() {
            assert(is_LPDDR4 && "Error: This functionality is not"
                    "currently supported for the selected DRAM standard.");

            if(num_grouped_SAs > 1)
                assert(num_weak_rows == 0 && "Error: Weak row remapping is not yet implemented to work with SA grouping");

            vector<int> addr_vec(int(T::Level::MAX), 0);


            for(int r = 0; r < spec->org_entry.count[int(T::Level::Rank)]; r++){
                addr_vec[int(T::Level::Rank)] = r;

                for(int b = 0; b < spec->org_entry.count[int(T::Level::Bank)]; b++) {
                    addr_vec[int(T::Level::Bank)] = b;

                    for(int i = 0; i < num_SAs; i++) {
                        addr_vec[int(T::Level::Row)] = i*SA_size;
                        int offset = calc_entries_offset(addr_vec);
                        for(int j = 0; j < num_weak_rows; j++) {
                            CROWEntry* entry = entries + offset + j;
                            entry->is_to_remap_weak_row = true;
                            entry->valid = true;
                        }
                    }
                }
            }
        }
 
};


} // namespace ramulator

#endif //__H_CROW_TABLE_H_