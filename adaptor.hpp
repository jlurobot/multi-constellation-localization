#pragma once

#include "nanoflann.hpp"

#include <vector>
#include <span>

// ===== This example shows how to use nanoflann with these types of containers: =======
//typedef std::vector<std::vector<double> > my_vector_of_vectors_t;
//typedef std::vector<Eigen::VectorXd> my_vector_of_vectors_t;   // This requires #include <Eigen/Dense>
// =====================================================================================


/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the storage.
  *  The i'th vector represents a point in the state space.
  *
  *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality for the points in the data set, allowing more compiler optimizations.
  *  \tparam num_t The type of the point coordinates (typically, double or float).
  *  \tparam Distance The distance metric to use: nanoflann::metric_L1, nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
  *  \tparam IndexType The type for indices in the KD-tree index (typically, size_t of int)
  */
template <typename point_type>
struct array_adaptor
{
    using num_t = decltype(point_type::x);
	typedef array_adaptor<point_type> self_t;
	typedef typename nanoflann::metric_L2::template traits<num_t,self_t>::distance_t metric_t;
	typedef nanoflann::KDTreeSingleIndexAdaptor<metric_t,self_t, 3 ,size_t>  index_t;

	index_t* index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.

	/// Constructor: takes a const ref to the vector of vectors object with the data points
	array_adaptor(const point_type* array, size_t count, const int leaf_max_size = 10): m_data(array), length(count)
	{
		index = new index_t(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size));
        index->buildIndex();
	}

	~array_adaptor() {
		delete index;
	}

	const point_type* m_data;
    size_t length;

	/** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
	  *  Note that this is a short-cut method for index->findNeighbors().
	  *  The user can also call index->... methods as desired.
	  * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
	  */
	template<typename query_point_type>
	inline void query(const query_point_type &query_point, const size_t num_closest, size_t *out_indices, num_t *out_distances_sq, const int nChecks_IGNORED = 10) const
	{   
        num_t val[3] = {query_point.x, query_point.y, query_point.z};
		nanoflann::KNNResultSet<num_t, size_t> resultSet(num_closest);
		resultSet.init(out_indices, out_distances_sq);
		index->findNeighbors(resultSet, val, nanoflann::SearchParams());
	}

	/** @name Interface expected by KDTreeSingleIndexAdaptor
	  * @{ */

	const self_t & derived() const {
		return *this;
	}
    
	self_t & derived()       {
		return *this;
	}

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const {
		return length;
	}

	// Returns the dim'th component of the idx'th point in the class:
	inline num_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        switch(dim) {
            case 0: return m_data[idx].x;
            case 1: return m_data[idx].y;
            case 2: return m_data[idx].z;
        }
        return 0;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX & /*bb*/) const {
		return false;
	}
}; // end of KDTreeVectorOfVectorsAdaptor