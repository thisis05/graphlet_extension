#include <cstddef>

#ifdef _WIN32
    using difference_type = std::ptrdiff_t;
#else
    #include <sys/types.h>
    using difference_type = ssize_t;
#endif

#include <iterator>

template <class Key, class Val>
struct JSIterator
{
  Key *key;
  Val *val;

  struct Pair
  {
    Key key;
    Val val;

    operator Key () const { return key; }
    bool operator < (const Pair &other) const { return key < other.key; }
    bool operator == (const Pair &other) const { return key == other.key; }
  };

  using value_type = Pair;
  using difference_type = difference_type; 
  using pointer = Pair*;
  using reference = Pair&;
  using iterator_category = std::random_access_iterator_tag;

  struct ReferenceProxy
  {
    Key *key;
    Val *val;
    
    ReferenceProxy& operator =(const value_type& v)
    {
      *key = v.key;
      *val = v.val;
      return *this;
    }

    ReferenceProxy& operator =(const ReferenceProxy& other)
    {
      *key = *other.key;
      *val = *other.val;
      return *this;
    }

    operator value_type() const
    {
      return {*key, *val};
    }

    operator Key () const { return *key; }

    bool operator < (const ReferenceProxy& other) const { return *key < *other.key; }
    bool operator < (const value_type& v) const { return *key < v.key; }

    friend void swap(ReferenceProxy a, ReferenceProxy b)
    {
      std::swap(*a.key, *b.key);
      std::swap(*a.val, *b.val);
    }
  };

  inline JSIterator& operator ++()
  {
    ++key;
    ++val;
    return *this;
  }

  inline JSIterator& operator --()
  {
    --key;
    --val;
    return *this;
  }

  inline JSIterator& operator +=(difference_type n)
  {
    key += n;
    val += n;
    return *this;
  }

  inline bool operator != (const JSIterator& other) const
  {
    return key != other.key;
  }

  inline bool operator == (const JSIterator& other) const
  {
    return key == other.key;
  }

  inline bool operator < (const JSIterator& other) const
  {
    return key < other.key;
  }

  inline bool operator > (const JSIterator& other) const
  {
    return key > other.key;
  }

  inline bool operator >= (const JSIterator& other) const
  {
    return key >= other.key;
  }

  JSIterator operator +(difference_type n) const
  {
    return {key + n, val + n};
  }

  JSIterator operator -(difference_type n) const
  {
    return {key - n, val - n};
  }

  difference_type operator -(const JSIterator& other) const
  {
    return key - other.key;
  }

  ReferenceProxy operator *() const { return {key, val}; }

  ReferenceProxy operator [](difference_type n) const { return {key + n, val + n}; }
};