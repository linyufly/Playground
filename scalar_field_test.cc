// Author: Mingcheng Chen (linyufly@gmail.com)

#include "util.h"

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

#include <queue>

namespace {

// Input files.
const char *kScalarFieldFile = "smoothed_scalar.vtk";
const char *kBinaryFieldFile = "binary_volume.vtk";
const char *kBinaryFieldToDeflate = "inflated_binary_volume.vtk";
const char *kMaskFieldFile = "mask_volume.vtk";

// Output files.
const char *kThresholdOutputFile = "threshold.vtk";
const char *kWaterproofOutputFile = "waterproof.vtk";
const char *kComponentOutputFile = "component.vtk";
const char *kInflateOutputFile = "inflate.vtk";
const char *kDeflateOutputFile = "deflate.vtk";
const char *kMaskOutputFile = "masked.vtk";

const double kThreshold = 4.5;  // 4.0, 4.5, 5.0

const int kConnectivity = 6;
const int kNumberOfInflations = 15;
const int kNumberOfDeflations = 15;

const int dire[6][3] = {
    {1, 0, 0}, {-1, 0, 0},
    {0, 1, 0}, {0, -1, 0},
    {0, 0, 1}, {0, 0, -1}};

int code(int *index, int *dimensions) {
  return (index[2] * dimensions[1] + index[1]) * dimensions[0] + index[0];
}

void decode(int code, int *dimensions, int *index) {
  index[0] = code % dimensions[0];
  code /= dimensions[0];
  index[1] = code % dimensions[1];
  index[2] = code / dimensions[1];
}

bool outside(int *index, int *dimensions) {
  for (int c = 0; c < 3; c++) {
    if (index[c] < 0 || index[c] >= dimensions[c]) {
      return true;
    }
  }

  return false;
}

int flood_fill(int s_x, int s_y, int s_z, vtkStructuredPoints *field,
                int color, int ***mark) {
  int dimensions[3];
  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int start_index[] = {s_x, s_y, s_z};
  std::queue<int> queue;
  queue.push(code(start_index, dimensions));

  mark[s_x][s_y][s_z] = color;

  int num_elements = 0;

  while (!queue.empty()) {
    num_elements++;

    int curr_code = queue.front();
    queue.pop();

    int curr_index[3];

    decode(curr_code, dimensions, curr_index);

    for (int d = 0; d < kConnectivity; d++) {
      int x = curr_index[0] + dire[d][0];
      int y = curr_index[1] + dire[d][1];
      int z = curr_index[2] + dire[d][2];

      int next_index[] = {x, y, z};
      int next_code = code(next_index, dimensions);

      if (outside(next_index, dimensions)) {
        continue;
      }

      if (mark[x][y][z]) {
        continue;
      }

      if (field->GetPointData()->GetScalars()->GetTuple1(next_code)
          <= kThreshold) {
        continue;
      }

      mark[x][y][z] = color;
      queue.push(next_code);
    }
  }

  return num_elements;
}

void inflate(int *dimensions, int num_inflations, bool ***inside) {
  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int ***steps = create_3d_array<int>(nx, ny, nz);

  std::queue<int> queue;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        if (!inside[x][y][z]) {
          steps[x][y][z] = -1;
        } else {
          steps[x][y][z] = 0;
          int curr_index[] = {x, y, z};
          int curr_code = code(curr_index, dimensions);
          queue.push(curr_code);
        }
      }
    }
  }

  /// DEBUG ///
  printf("queue.size() = %d\n", static_cast<int>(queue.size()));

  while (!queue.empty()) {
    int curr_code = queue.front();
    queue.pop();

    int curr_index[3];
    decode(curr_code, dimensions, curr_index);

    for (int d = 0; d < kConnectivity; d++) {
      int x = curr_index[0] + dire[d][0];
      int y = curr_index[1] + dire[d][1];
      int z = curr_index[2] + dire[d][2];

      int next_index[] = {x, y, z};

      if (outside(next_index, dimensions)) {
        continue;
      }

      if (steps[x][y][z] != -1) {
        continue;
      }

      /// DEBUG ///
      // printf("x, y, z = %d, %d, %d\n", x, y, z);

      steps[x][y][z] = steps[curr_index[0]][curr_index[1]][curr_index[2]] + 1;
      int next_code = code(next_index, dimensions);

      queue.push(next_code);
    }
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        if (steps[x][y][z] <= num_inflations) {
          inside[x][y][z] = true;
        }

        if (steps[x][y][z] == -1) {
          printf("How could this happen? %d, %d, %d\n", x, y, z);
          exit(0);
        }
      }
    }
  }

  delete_3d_array(steps);
}

}

void threshold_test() {
  printf("threshold_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  for (int i = 0;
       i < field->GetPointData()->GetScalars()->GetNumberOfTuples();
       i++) {
    double value = field->GetPointData()->GetScalars()->GetTuple1(i);
    if (value <= kThreshold) {
      value = 0.0;
    } else {
      value = 1.0;
    }

    field->GetPointData()->GetScalars()->SetTuple1(i, value);
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kThresholdOutputFile);
  writer->SetInputData(field);
  writer->Write();

  printf("} threshold_test\n");
}

void waterproof_test() {
  printf("waterproof_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];

  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  bool ***visited = create_3d_array<bool>(nx, ny, nz);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        visited[x][y][z] = false;
      }
    }
  }

  visited[0][0][0] = true;
  std::queue<int> queue;
  queue.push(0);

  while (!queue.empty()) {
    int curr = queue.front();
    queue.pop();

    int curr_index[3];
    decode(curr, dimensions, curr_index);

    for (int d = 0; d < kConnectivity; d++) {
      int x = curr_index[0] + dire[d][0];
      int y = curr_index[1] + dire[d][1];
      int z = curr_index[2] + dire[d][2];

      int next_index[] = {x, y, z};
      if (outside(next_index, dimensions)) {
        continue;
      }

      if (visited[x][y][z]) {
        continue;
      }

      int next_code = code(next_index, dimensions);
      double value = field->GetPointData()->GetScalars()->GetTuple1(next_code);

      if (value > kThreshold) {
        continue;
      }

      visited[x][y][z] = true;
      queue.push(next_code);
    }
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);
        if (visited[x][y][z]) {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 0.0);
        } else {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 1.0);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kWaterproofOutputFile);
  writer->SetInputData(field);
  writer->Write();

  delete_3d_array(visited);

  printf("} waterproof_test\n");
}

void component_test() {
  printf("component_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];

  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int ***mark = create_3d_array<int>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        mark[x][y][z] = 0;
      }
    }
  }

  int num_colors = 0;
  int best_color = -1, largest = 0;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        if (!mark[x][y][z]) {
          int curr_index[] = {x, y, z};
          int curr_code = code(curr_index, dimensions);

          if (field->GetPointData()->GetScalars()->GetTuple1(curr_code)
              > kThreshold) {
            num_colors++;
            int num_elements = flood_fill(x, y, z, field, num_colors, mark);

            if (num_elements > largest) {
              largest = num_elements;
              best_color = num_colors;
            }

            printf("color %d: num_elements = %d (best color: %d, largest: %d)\n",
                   num_colors, num_elements, best_color, largest);
          }
        }
      }
    }
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        if (mark[x][y][z] != best_color) {
          mark[x][y][z] = 0;
        } else {
          mark[x][y][z] = 1;
        }

        field->GetPointData()->GetScalars()->SetTuple1(
            curr_code, static_cast<double>(mark[x][y][z]));
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kComponentOutputFile);
  writer->SetInputData(field);
  writer->Write();

  delete_3d_array(mark);

  printf("} component_test\n");
}

void inflate_test() {
  printf("inflate_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kBinaryFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];

  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  bool ***inside = create_3d_array<bool>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        double value = field->GetPointData()->GetScalars()
                                            ->GetTuple1(curr_code);

        if (value > 0.5) {
          inside[x][y][z] = true;
        } else {
          inside[x][y][z] = false;
        }
      }
    }
  }

  inflate(dimensions, kNumberOfInflations, inside);

  bool ***visited = create_3d_array<bool>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        visited[x][y][z] = false;
      }
    }
  }

  printf("inside[0][0][0] = %d\n", static_cast<int>(inside[0][0][0]));

  std::queue<int> queue;
  queue.push(0);
  visited[0][0][0] = true;

  while (!queue.empty()) {
    int curr_code = queue.front();
    queue.pop();

    int curr_index[3];
    decode(curr_code, dimensions, curr_index);

    for (int d = 0; d < kConnectivity; d++) {
      int x = curr_index[0] + dire[d][0];
      int y = curr_index[1] + dire[d][1];
      int z = curr_index[2] + dire[d][2];

      int next_index[] = {x, y, z};

      if (outside(next_index, dimensions)) {
        continue;
      }

      if (visited[x][y][z] || inside[x][y][z]) {
        continue;
      }

      int next_code = code(next_index, dimensions);
      visited[x][y][z] = true;
      queue.push(next_code);
    }
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        if (visited[x][y][z]) {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 0.0);
        } else {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 1.0);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kInflateOutputFile);
  writer->SetInputData(field);
  writer->Write();

  delete_3d_array(visited);
  delete_3d_array(inside);

  printf("} inflate_test\n");
}

void deflate_test() {
  printf("deflate_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kBinaryFieldToDeflate);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];

  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  bool ***outside = create_3d_array<bool>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        double value = field->GetPointData()->GetScalars()
                                            ->GetTuple1(curr_code);

        if (value > 0.5) {
          outside[x][y][z] = false;
        } else {
          outside[x][y][z] = true;
        }
      }
    }
  }

  inflate(dimensions, kNumberOfDeflations, outside);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        if (outside[x][y][z]) {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 0.0);
        } else {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 1.0);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kDeflateOutputFile);
  writer->SetInputData(field);
  writer->Write();

  delete_3d_array(outside);

  printf("} deflate_test\n");
}

void mask_test() {
  printf("mask_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  reader->SetFileName(kMaskFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> mask =
      vtkSmartPointer<vtkStructuredPoints>::New();
  mask->DeepCopy(reader->GetOutput());

  int dimensions[3];

  field->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);

        double mask_value = mask->GetPointData()->GetScalars()
                                                ->GetTuple1(curr_code);

        if (mask_value < 0.5) {
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, 0.0);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kMaskOutputFile);
  writer->SetInputData(field);
  writer->Write();

  printf("} mask_test\n");
}

int main() {
  // threshold_test();
  // waterproof_test();
  // component_test();
  // inflate_test();
  // deflate_test();
  mask_test();

  return 0;
}

