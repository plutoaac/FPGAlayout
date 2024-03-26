
#pragma once

template <class T, class Compare>
void quick_sort(T *data, std::size_t size, Compare compare) {
  if (size < 1) return;
  if (size < (1 << 12)) {
    std::sort(data, data + size, compare);
    return;
  }

  size_t mid = std::hash<size_t>{}(size);
  mid ^= std::hash<void *>{}(static_cast<void *>(data));
  mid %= size;
  std::swap(data[0], data[mid]);
  T pivot = data[0];
  size_t left = 0, right = size - 1;
  while (left < right) {
    while (left < right && !compare(data[right], pivot)) right--;
    if (left < right) data[left++] = data[right];
    while (left < right && compare(data[left], pivot)) left++;
    if (left < right) data[right--] = data[left];
  }
  data[left] = pivot;

#pragma omp parallel sections
  {
#pragma omp section
    quick_sort(data, left, compare);

#pragma omp section
    quick_sort(data + left + 1, size - left - 1, compare);
  }
}

template <class T>
void quick_sort(T *data, int size) {
  quick_sort(data, size, std::less<T>());
}