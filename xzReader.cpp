#include <iostream>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include <fcntl.h>
#include <vector>
#include <lzma.h>
#include <string.h>
#include <Rcpp.h>


typedef uint16_t vbin_t;
static const long double VBIN_T_MAX = (pow(2,8*sizeof(vbin_t)) - 1);

static const bool xz_mode = true;
static const size_t XZ_BUFFER_SIZE =  2 * 1024 * 1024;

static const int mac_dec = 6;
static const double mac_tol = 2.00 * std::pow(0.1, (double) mac_dec);

bool treat_as_int(const double& m){
	return ( std::abs(m - round(m)) < mac_tol );	
}

double unflip(const int& f1, const int& f2){
        if( f1 != f2 ){
                return -1.0;
        }else{
                return 1.0;
        }
}

unsigned int pack_dp(const double& x, const double& m0, const double& m1){
	if( m0 > m1 ){
		return pack_dp(x, m1, m0);
	}else{
		if( m0 <= 0.0 ){
			return (unsigned int) 0;
		}else if( 2.0 * m0 < VBIN_T_MAX && treat_as_int(m0) && treat_as_int(m1) ){
			return (unsigned int) round(x);
		}else{
			return (unsigned int) round(VBIN_T_MAX * x/ceil((2.0*m0)));
		}
	}
}

double unpack_dp(const double& x, const double& m0, const double& m1){
	if( m0 > m1 ){
		return unpack_dp(x, m1, m0);
	}else{
		if( m0 <= 0.0 ){
			return 0.0;
		}else if( 2.0 * m0 < VBIN_T_MAX && treat_as_int(m0) && treat_as_int(m1)){
			return x;
		}else{
			return x*ceil(2.0*m0)/((double) VBIN_T_MAX);
		}
	}
}


class call_tracker
{
        private:
                std::vector<int> most_recent;
                int n;
        public:
                call_tracker(const int& size) : n(size) {};
                call_tracker() {};
                void set_size(const int& n_recent){
                        n = n_recent;
                }

                int push_new(const int& x){
                        most_recent.push_back(x);
                        if( most_recent.size() < n ){
                                // Don't remove anything from cache.
                                return -1;
                        }
                        int out = most_recent.front();
                        most_recent.erase(most_recent.begin());
                        return out;
                }

                int check_new( const int& x ){
                        if( most_recent.size() < n ){
                                most_recent.push_back(x);
                                return -1;
                        }
                        int x_i = -1;
                        for(int i = n; i >= 0; i-- ){
                                if( most_recent[i] == x ){
                                        x_i = i;
                                        break;
                                }
                        }
                        if( x_i > 0 ){
                                std::swap(most_recent[most_recent.size()-1], most_recent[x_i]);
                                return -1;
                        }else{
                                most_recent.push_back(x);
                                int out = most_recent.front();
                                most_recent.erase(most_recent.begin());
                                // Return the value that we want to delete from cache.
                                return out;
                        }
                }
};


// Code below for processing xz files is partly based on the
// following GitHub repository:
//     https://github.com/CTSRD-CHERI/cheritrace.git
// The license issued by the original authors
// (Alfredo Mazzinghi and David T. Chisnall) for this code is here:
//     http://www.beri-open-systems.org/legal/license-1-0.txt

class xzReader
{
        protected:
                struct offsets
                {
                        off_t c_start; size_t c_size;
                        off_t u_start; size_t u_size;
                };
                size_t read_compressed(void *buffer, off_t start, size_t n_bytes){
                        if ( start < 0 || start > compressed_file_size ){
                                return 0;
                        }
                        if (start + n_bytes > compressed_file_size){
                                n_bytes = compressed_file_size - start;
                        }
                        size_t n_completed = 0;
                        while (n_bytes > 0)
                        {
                                ssize_t n_bytes_i = pread(compressed_file, buffer, n_bytes, start);
                                if (n_bytes_i < 0) break;
                                n_completed += n_bytes_i;
                                start += n_bytes_i;
                                n_bytes -= n_bytes_i;
                                buffer = (void*)((char*)buffer + n_bytes_i);
                        }
                        return n_completed;
                }

                call_tracker cache_czar;
                std::vector<std::unique_ptr<uint8_t>> cache;

                int compressed_file;

                std::vector<offsets> block_offsets;
                size_t compressed_file_size = 0;

                lzma_stream_flags stream_flags;

                size_t read(void *buffer, off_t start, size_t length)
                {
                        int block_idx = get_block_from_offset(start);
                        if (block_idx < 0)
                        {
                                return 0;
                        }
                        size_t copied = 0;
                        while (length > 0)
                        {
                                if (block_idx >= (int)block_offsets.size())
                                {
                                        break;
                                }
                                if( cache_block(block_idx) < 0 ){
                                        Rcpp::Rcout << "ERROR: Cannot read block " << block_idx << "\n";
                                        return 0;
                                }
                                const auto& data = cache[block_idx];
                                auto &b = block_offsets[block_idx++];
                                size_t copy_start = start - b.u_start;
                                size_t copy_length = b.u_size - copy_start;
                                copy_length = std::min(copy_length, length);
                                memcpy(buffer, data.get()+copy_start, copy_length);
                                copied += copy_length;
                                start += copy_length;
                                length -= copy_length;
                                buffer = (void*)((char*)buffer + copy_length);
                        }
                        return copied;
                }

                int get_block_from_offset(off_t off)
                {
                        int i = 0;
                        for(const offsets& bi : block_offsets){
                                if( (bi.u_start <= off) && (bi.u_start + bi.u_size > off) ){
                                        return i;
                                }
                                i++;
                        }
                        Rcpp::Rcout << "Warning: Requested data outside valid range.\n";
                        return -1;
                }

                int cache_block(const int& i)
                {
                        if( cache[i] != nullptr ){
                                return 1;
                        }else{
                                int prune_cache = cache_czar.push_new(i);
                                if( prune_cache > 0 ){
                                        cache[prune_cache] = nullptr;
                                }
                        }

                        offsets& b = block_offsets[i];

                        std::unique_ptr<uint8_t> buffer(new uint8_t[b.c_size]);

                        read_compressed((void*)buffer.get(), b.c_start, b.c_size);
                        lzma_block block;
                        lzma_filter filters[LZMA_FILTERS_MAX + 1];

                        filters[0].id = LZMA_VLI_UNKNOWN;
                        block.filters = filters;
                        block.version = 1;
                        block.check = stream_flags.check;
                        block.header_size = lzma_block_header_size_decode(*buffer);

                        if ( lzma_block_header_decode(&block, nullptr, buffer.get()) != LZMA_OK )
                        {
                                return -1;
                        }

                        cache[i] = std::unique_ptr<uint8_t>(new uint8_t[b.u_size]);
                        size_t in_pos = block.header_size;
                        size_t out_pos = 0;

                        if ( lzma_block_buffer_decode(&block, nullptr, buffer.get(),&in_pos, b.c_size, cache[i].get(), &out_pos, b.u_size) != LZMA_OK)
                        {
                                return -1;
                        }

                        return 1;
                }

        public:

                void open(std::string filepath)
                {

                        compressed_file = ::open(filepath.c_str(), O_RDONLY);
                        compressed_file_size = lseek(compressed_file, 0, SEEK_END);

                        uint8_t footer_buffer[12];
                        read_compressed((void*)footer_buffer, compressed_file_size-12, 12);

                        if ( lzma_stream_footer_decode(&stream_flags, footer_buffer) != LZMA_OK ){
                                // cerr << "ERROR: Cannot read xz footer.\n";
                                close(compressed_file);
                                abort();
                        }

                        std::unique_ptr<uint8_t> index_buffer(new uint8_t[stream_flags.backward_size]);
                        read_compressed((void*)index_buffer.get(),
                                                                  compressed_file_size - stream_flags.backward_size - 12,
                                                                  stream_flags.backward_size);
                        lzma_index *idx;
                        uint64_t mem = UINT64_MAX;
                        size_t pos = 0;

                        if ( lzma_index_buffer_decode(&idx, &mem, nullptr,index_buffer.get(), &pos, stream_flags.backward_size) != LZMA_OK ){
                                // cerr << "ERROR: Cannot read xz index.\n";
                                close(compressed_file);
                                abort();
                        }
                        lzma_index_iter iter;
                        lzma_index_iter_init(&iter, idx);
                        int bl = 0;
                        while (!lzma_index_iter_next(&iter, LZMA_INDEX_ITER_ANY))
                        {
                                offsets block;
                                block.c_start = iter.block.compressed_file_offset;
                                block.c_size = iter.block.total_size;
                                block.u_start = iter.block.uncompressed_file_offset;
                                block.u_size = iter.block.uncompressed_size;
                                block_offsets.push_back(block);
                                bl++;
                        }
                        lzma_index_end(idx, nullptr);

                        cache_czar.set_size(10);
                        cache.resize(block_offsets.size());
                }
                xzReader(std::string filepath){
                        open(filepath);
                }
                xzReader() {};
                ~xzReader() { cache.clear(); };
                Rcpp::NumericVector getData(int64_t st, int64_t n_vals){

                        int64_t bytes_read = 1;
                        int64_t start, end, size;

                        start = st;
                        size = -1; end = -1;

                        int64_t nval_buffer;
                        end = start + n_vals*sizeof(vbin_t);
                        void* buffer = malloc(XZ_BUFFER_SIZE);
                        Rcpp::NumericVector out_vec(n_vals);

                        int64_t ii = 0;

                        while( start <= end && ii < n_vals){

                                int64_t nbytes = XZ_BUFFER_SIZE;

                                if (end >= 0){
                                        nbytes = (end - start > XZ_BUFFER_SIZE)? XZ_BUFFER_SIZE:(end - start);
                                }
                                bytes_read = read(buffer, (off_t) start, (size_t) nbytes);

                                if (bytes_read == 0) break;
                                if (bytes_read < 0){
                                        Rcpp::Rcout << "ERROR reading from VCOV file.\n";
                                        break;
                                }
                                start += bytes_read;

                                nval_buffer = (int)bytes_read/sizeof(vbin_t);

                                for( int i = 0; i < nval_buffer; i++ ){
                                        if( ii >= n_vals ){
                                                break;
                                        }else{
                                                vbin_t value = *reinterpret_cast<vbin_t*>((char*)buffer+sizeof(vbin_t)*i);
                                                out_vec[ii] = (double) value;
                                                ii++;
                                        }
                                }
                        }

                        free(buffer);
                        buffer = NULL;

                        return out_vec;
                }
};

using namespace Rcpp;

NumericMatrix buildMatrixC(const NumericVector& s, const NumericVector& n, const NumericVector& m, const NumericVector& x){
        int64_t sz = s.size();
        int64_t sz_x = x.size();
        int64_t s0 = s[0];
        NumericMatrix out(sz, sz);
        for(int64_t i = 0; i < sz; i++ ){
                if( m[i] > 0 ){
			int64_t offset_ij = s[i] - s0;
                        if( offset_ij >= sz_x ){
                                Rf_error("ERROR 1");
                        }
			out(i,i) = unpack_dp(x[offset_ij], m[i], m[i]);
			
                        int64_t n_i = i + n[i];
                        n_i = n_i < sz ? n_i : sz;
                        if( i < sz - 1 ){
                                for(int64_t j = i + 1; j < n_i; j++){
					offset_ij++;
                                        if( offset_ij >= sz_x ){
                                                Rf_error("ERROR 2");
                                        }
                                        out(j,i) = unpack_dp(x[offset_ij], m[i], m[j]);
					out(i,j) = out(j,i);
                                }
                        }
                }
        }
        return out;
}

NumericVector buildVectorC(const int64_t& target, const NumericVector& s, const NumericVector& n, const NumericVector& m, const NumericVector& x){
        int64_t sz = s.size();
        int64_t sz_x = x.size();
        int64_t s0 = s[0];
	double m_t = m[target];
	int64_t s_t = s[target];
	int64_t n_t = n[target];

        NumericVector out(sz);
		
	if( m_t <= 0 ){
		return out;
	}

	for(int64_t i = 0; i < target; i++ ){
		if( m[i] > 0 && n[i] > target - i ){
			int64_t offset_i = s[i] - s0 + (target - i);
			if( offset_i >= sz_x ){
				Rf_error("ERROR: requested LD outside window.");
			}
			out(i) = unpack_dp(x[offset_i], m[i], m_t);
		}
	}
	int64_t offset_i = s_t;
	for(int64_t i = target; i < sz; i++ ){
		if( offset_i >= sz_x || i - target >= n_t ){
			Rf_error("ERROR: requested LD outside window.");
		}
		out(i) = unpack_dp(x[offset_i], m[i], m_t);
		offset_i++;
	}
        return out;
}


template <typename T>
bool all_lt( const std::vector<T>& ii, const std::vector<T>& nn ){
        for(int i = 0; i < ii.size(); ++i){
                if( ii[i] >= nn[i] ){
                        return false;
                }
        }
        return true;
}

template <typename T>
void seq_to(std::vector<T>& ii, const std::vector<T>& n_var, const std::vector<std::vector<T>>& pos, const T& target){
        for( int s = 0; s < ii.size(); s++ ){
                while( pos[s][ii[s]] < target && ii[s] < n_var[s] ){
                        ii[s]++;
                }
        }
}

std::vector<std::vector<int64_t>> mergeIntersect(const std::vector<std::vector<int64_t>>& pos, const std::vector<StringVector>& ref, const std::vector<StringVector>& alt){
        int64_t N = pos.size();
        std::vector<std::vector<int64_t>> out(N);
        std::vector<int64_t> n_var(N);
        std::vector<int64_t> ii(N, 0);
        for( int i = 0 ; i < N; i++){
                n_var[i] = pos[i].size();
        }
        while( all_lt(ii, n_var) ){
                int mp = 0;
                bool skip = false;

                for( int s = 0; s < N; s++ ){
                        // Rcout << ii[s] << ", ";
                        if( pos[s][ii[s]] > mp ){
                                mp = pos[s][ii[s]];
                        }
                }
                // Rcout << "\n";
                for( int s = 0; s < N; s++ ){
                        while( pos[s][ii[s]] < mp && ii[s] < n_var[s] ){
                                ii[s]++;
                        }

                        if( ii[s] >= n_var[s] ){
                                break;
                        }else if( pos[s][ii[s]] != mp ){
                                skip = true;
                        }else if( ii[s] + 1 < n_var[s] ){
                                if( pos[s][ii[s] + 1] == mp ){
                                        skip = true;
                                }
                        }
                }
                if( !all_lt(ii, n_var) ) break;
                if( skip ){
                        for( int s = 0; s < N; s++ ){
                                while( pos[s][ii[s]] <= mp && ii[s] < n_var[s] ){
                                        ii[s]++;
                                }
                        }
                        continue;
                }

                for( int s = 1; s < N; s++ ){
                        if( !( ref[s][ii[s]] == ref[0][ii[0]] && alt[s][ii[s]] == alt[0][ii[0]] ) ){
                                ii[s]++;
                                skip = true;
                        }
                }
                if( skip ) continue;

                for( int s = 0; s < N; s++ ){
                        out[s].push_back(ii[s] + 1);
                        ii[s]++;
                }
        }

        return out;
}

NumericMatrix flipMatrix(NumericMatrix A, const NumericVector& w_flip, const NumericVector& w_good){
        int n = A.cols();
        for(const auto& i : w_good){
                if( i < 1 || i > n ){
                        Rf_error("ERROR");
                }
        }
        for(const auto& i : w_flip){
                if( i < 1 || i > n ){
                        Rf_error("ERROR");
                }
                for(const auto& j : w_good){
                        A(i-1,j-1) *= (-1);
                        A(j-1,i-1) *= (-1);
                }
        }
        return(A);
};

//void finalize_xzReader( xzReader* pt ){
//        pt->finalize();
//        return;
//};


RCPP_EXPOSED_CLASS(xzReader)
RCPP_MODULE(mod_test) {

        function("buildMatrixC", &buildMatrixC);
        function("flipMatrix", &flipMatrix);
        function("mergeIntersect", &mergeIntersect);

    class_<xzReader>("xzReader")

        .constructor<std::string>("Open an xz file for random access of uncompressed data.")

        .method("open", &xzReader::open, "Open an xz file.")
        .method("getData", &xzReader::getData, "Random access to uncompressed data from xz file.")
        // .finalizer( &finalize_xzReader )
        ;

}
