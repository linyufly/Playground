// Author: Mingcheng Chen (linyufly@gmail.com)

#include "util.h"

#include <cmath>

#include <algorithm>
#include <map>
#include <queue>
#include <vector>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

namespace {

// Input files.
// const char *kScalarFieldFile = "smoothed_scalar.vtk";
// const char *kScalarFieldFile = "smoothed_scalar_5.vtk";
// const char *kScalarFieldFile = "/home/linyufly/Data/P96_bFTLE.vtk";
// const char *kScalarFieldFile = "/home/linyufly/GitHub/GaussianSmoothing3D/gaussian_smoothed.vtk";
// const char *kScalarFieldFile = "P96_bFTLE/FTLE96ixbwd3D_highres_noT.1.vtk";
// const char *kScalarFieldFile = "/home/linyufly/GitLab/RidgeExtraction/data/convective_half_ftle.vtk";
const char *kScalarFieldFile = "/home/linyufly/GitLab/RidgeExtraction/data/convective_half_label_dig_corrected.vtk";
// const char *kScalarFieldFile = "P96_bFTLE/FTLE96ixbwd3D_highres_noT.1.vtk";
// const char *kScalarFieldFile = "P96_bFTLE/FTLE96ixbwd3D_highres_noT.2.vtk";
// const char *kScalarFieldFile = "P96_bFTLE/FTLE96ixbwd3D_highres_noT.24.vtk";
// const char *kBinaryFieldFile = "binary_volume.vtk";
// const char *kBinaryFieldFile = "threshold_volume_raw_5_2.vtk";
// const char *kBinaryFieldFile = "largest_component.vtk";
const char *kBinaryFieldFile = "component.vtk";
// const char *kBinaryFieldToDeflate = "inflated_binary_volume.vtk";
// const char *kBinaryFieldToDeflate = "inflated_binary_volume_raw_5_2.vtk";
const char *kBinaryFieldToDeflate = "inflate.vtk";
// const char *kMaskFieldFile = "mask_volume.vtk";
// const char *kMaskFieldFile = "deflated_binary_volume_raw_5_2.vtk";
const char *kMaskFieldFile = "deflate.vtk";
// const char *kMaskFieldFile = "inflated_binary_volume.vtk";
// const char *kMaskFieldFile = "deflate.vtk";
// const char *kCoarseLabelFile = "coarse_labelling.vtk";
const char *kCoarseLabelFile = "/home/linyufly/GitHub/Watershed/itk_watershed_image_filter_result.vtk";
// const char *kFineLabelFile = "fine_labelling.vtk";
const char *kFineLabelFile = "/home/linyufly/GitHub/Watershed/itk_watershed_image_filter_result.vtk";
const char *kScalarFieldToCrop = "/home/linyufly/Data/convective_half.vtk";
const char *kScalarFieldToSymmetricalize = "cropped_convective_half.vtk";
// const char *kLabelFieldToModify = "convective_labelling.vtk";
const char *kLabelFieldToModify = "modify.vtk";
const char *kFineLabelFieldFile = "convective_half_fine_labelling.vtk";
const char *kCoarseLabelFieldFile = "convective_half_coarse_labelling.vtk";

// Output files.
const char *kThresholdOutputFile = "threshold.vtk";
const char *kWaterproofOutputFile = "waterproof.vtk";
const char *kComponentOutputFile = "component.vtk";
const char *kInflateOutputFile = "inflate.vtk";
const char *kDeflateOutputFile = "deflate.vtk";
const char *kMaskOutputFile = "masked.vtk";
const char *kBoundaryOutputFile = "boundary.vtk";
const char *kCropOutputFile = "crop.vtk";
const char *kSymmetricalizeOutputFile = "symmetricalize.vtk";
const char *kModifyOutputFile = "modify.vtk";
const char *kDigOutputFile = "dig.vtk";
const char *kMirrorOutputFile = "mirror.vtk";

const double kThreshold = 5.8;  // 4.0, 4.5, 5.0, 5.5, 6.0
                                // 5.8 for 0
                                // 5.8 for 1
                                // 5.7 for 2
                                // 5.4 for 3
                                // 5.4 for 4
                                // 5.4 for 5
                                // 5.4 for 6
                                // 5.3 for 7
                                // 5.3 for 8
                                // 5.3 for 9
                                // 5.2 for 10
                                // 5.2 for 11
                                // 5.2 for 12
                                // 5.3 for 13
                                // 5.3 for 14
                                // 5.4 for 15
                                // 5.4 for 16
                                // 5.4 for 17
                                // 5.4 for 18
                                // 5.4 for 19
                                // 5.4 for 20
                                // 5.4 for 21
                                // 5.5 for 22
                                // 5.7 for 23
                                // 5.7 for 24
const double kCropDistance = 1.6;

const int kConnectivity = 6;
const int kNumberOfInflations = 10;
const int kNumberOfDeflations = 5;

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

double square(double value) {
  return value * value;
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
          field->GetPointData()->GetScalars()->SetTuple1(curr_code, -1.0);
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

void boundary_test() {
  printf("boundary_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kFineLabelFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> fine =
      vtkSmartPointer<vtkStructuredPoints>::New();
  fine->DeepCopy(reader->GetOutput());

  reader->SetFileName(kCoarseLabelFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> coarse =
      vtkSmartPointer<vtkStructuredPoints>::New();
  coarse->DeepCopy(reader->GetOutput());

  int dimensions[3];

  fine->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  std::map<int, int> label_map;
  int num_labels = 0;

  int ***labels = create_3d_array<int>(nx, ny, nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);
        double curr_label =
            coarse->GetPointData()->GetScalars()->GetTuple1(curr_code);

        int value = static_cast<int>(curr_label + 0.5);

        if (label_map.find(value) == label_map.end()) {
          label_map[value] = num_labels++;
        }

        labels[x][y][z] = label_map[value];
      }
    }
  }

  bool detected = false;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        /// DEBUG ///
        // if (labels[x][y][z] != labels[0][0][0]) {
        //   continue;
        // }

        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);
        double label_0 =
            fine->GetPointData()->GetScalars()->GetTuple1(0);
        double curr_label =
            fine->GetPointData()->GetScalars()->GetTuple1(curr_code);

        if (curr_label != label_0) {
          labels[x][y][z] = num_labels;
          detected = true;
        }
      }
    }
  }

  if (detected) {
    num_labels++;
  }

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int curr_index[] = {x, y, z};
        int curr_code = code(curr_index, dimensions);
        coarse->GetPointData()->GetScalars()->SetTuple1(
            curr_code, static_cast<double>(labels[x][y][z]));
      }
    }
  }

  delete_3d_array(labels);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kBoundaryOutputFile);
  writer->SetInputData(coarse);
  writer->Write();

  printf("} boundary_test\n");
}

void crop_test() {
  printf("crop_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldToCrop);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];
  field->GetDimensions(dimensions);

  double spacing[3], origin[3];
  field->GetSpacing(spacing);
  field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;

        double coord[3];
        coord[0] = origin[0] + spacing[0] * x;
        coord[1] = origin[1] + spacing[1] * y;
        coord[2] = origin[2] + spacing[2] * z;

        if (sqrt(coord[0] * coord[0] + coord[1] * coord[1]) > kCropDistance) {
          field->GetPointData()->GetScalars()->SetTuple1(index, 0.0);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kCropOutputFile);
  writer->SetInputData(field);
  writer->Write();

  printf("} crop_test\n");
}

void symmetricalize_test() {
  printf("symmetricalize_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldToSymmetricalize);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];
  field->GetDimensions(dimensions);

  double spacing[3], origin[3];
  field->GetSpacing(spacing);
  field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  for (int z = 0; z < nz; z++) {
    printf("z = %d\n", z);

    std::map<double, double> sum_value;
    std::map<double, int> count;

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        int index = (z * ny + y) * nx + x;

        double coord[2];
        coord[0] = origin[0] + spacing[0] * x;
        coord[1] = origin[1] + spacing[1] * y;

        double dist = sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
        double value = field->GetPointData()->GetScalars()->GetTuple1(index);

        sum_value[dist] += value;
        count[dist]++;
      }
    }

    std::vector<double> distance_list;
    for (std::map<double, double>::iterator itr = sum_value.begin();
         itr != sum_value.end(); ++itr) {
      distance_list.push_back(itr->first);
    }

    // double tolerance = sqrt(spacing[0] * spacing[0] + spacing[1] * spacing[1]);
    double tolerance = spacing[0];

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        int index = (z * ny + y) * nx + x;

        double coord[2];
        coord[0] = origin[0] + spacing[0] * x;
        coord[1] = origin[1] + spacing[1] * y;

        double dist = sqrt(coord[0] * coord[0] + coord[1] * coord[1]);

        int lb = lower_bound(distance_list.begin(), distance_list.end(),
                             dist - tolerance) - distance_list.begin();
        int ub = upper_bound(distance_list.begin(), distance_list.end(),
                             dist + tolerance) - distance_list.begin();

        double sum_scalar = 0.0;
        int sum_count = 0;

        for (int p = lb; p < ub; p++) {
          double curr_dist = distance_list[p];
          if (sum_value.find(curr_dist) == sum_value.end()) {
            continue;
          }

          sum_scalar += sum_value[curr_dist];
          sum_count += count[curr_dist];
        }

        /// DEBUG ///
        if (sum_count == 0) {
          printf("dist = %lf\n", dist);
          printf("x, y, z = %d, %d, %d\n", x, y, z);
          printf("lb, ub = %d, %d\n", lb, ub);
          for (int p = lb; p < ub; p++) {
            printf(" %lf", distance_list[p]);
          }
          printf("\n");
          exit(0);
        }

        double average = sum_scalar / sum_count;

        field->GetPointData()->GetScalars()->SetTuple1(index, average);
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kSymmetricalizeOutputFile);
  writer->SetInputData(field);
  writer->Write();

  printf("} symmetricalize_test\n");
}

void modify_test() {
  printf("modify_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kLabelFieldToModify);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  field->DeepCopy(reader->GetOutput());

  int dimensions[3];
  field->GetDimensions(dimensions);

  double spacing[3], origin[3];
  field->GetSpacing(spacing);
  field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  double coord_1[] = {0.0, 1.0, 0.9};
  double coord_2[] = {0.0, 1.0, 2.0};

  int y_1 = static_cast<int>((coord_1[1] - origin[1]) / spacing[1]);
  int z_1 = static_cast<int>((coord_1[2] - origin[2]) / spacing[2]);

  int y_2 = static_cast<int>((coord_2[1] - origin[1]) / spacing[1]);
  int z_2 = static_cast<int>((coord_2[2] - origin[2]) / spacing[2]);

  int index_1 = (z_1 * ny + y_1) * nx;
  int index_2 = (z_2 * ny + y_2) * nx;

  double label_1 = field->GetPointData()->GetScalars()->GetTuple1(index_1);
  double label_2 = field->GetPointData()->GetScalars()->GetTuple1(index_2);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        if (field->GetPointData()->GetScalars()->GetTuple1(index) == label_1) {
          field->GetPointData()->GetScalars()->SetTuple1(index, label_2);
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetFileName(kModifyOutputFile);
  writer->SetInputData(field);
  writer->Write();
 

  printf("} modify_test\n");
}

void dig_test() {
  printf("dig_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kFineLabelFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> fine_label_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  fine_label_field->DeepCopy(reader->GetOutput());

  reader->SetFileName(kCoarseLabelFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> coarse_label_field =
      vtkSmartPointer<vtkStructuredPoints>::New();
  coarse_label_field->DeepCopy(reader->GetOutput());

  int dimensions[3];
  coarse_label_field->GetDimensions(dimensions);

  double origin[3], spacing[3];
  coarse_label_field->GetOrigin(origin);
  coarse_label_field->GetSpacing(spacing);

  double points[4][3] = {{0.0, 1.15, 1.66},
                         {0.0, 1.15, 1.56},
                         {0.0, 1.08, 1.48},
                         {0.0, 1.08, 1.42}};

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int point_indices[4][3];
  for (int p = 0; p < 4; p++) {
    for (int c = 0; c < 3; c++) {
      point_indices[p][c] = static_cast<int>((points[p][c] - origin[c]) / spacing[c]);
    }
  }

  double point_labels[4];
 
  int first_index;
  for (int p = 0; p < 4; p++) {
    int index = (point_indices[p][2] * ny + point_indices[p][1]) * nx + point_indices[p][0];

    if (!p) {
      first_index = index;
    }

    point_labels[p] = fine_label_field->GetPointData()->GetScalars()->GetTuple1(index);
  }

  int num_labels = 0;
  std::map<double, int> label_map;

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;

        double label = coarse_label_field->GetPointData()->GetScalars()->GetTuple1(index);

        if (label_map.find(label) == label_map.end()) {
          label_map[label] = num_labels++;
        }

        coarse_label_field->GetPointData()->GetScalars()->SetTuple1(index, label_map[label]);
      }
    }
  }

  double original_label = coarse_label_field->GetPointData()->GetScalars()->GetTuple1(first_index);

  for (int z = 0; z < nz; z++) {
    double min_dist = -1.0;
    double max_dist = -1.0;

    double curr_z = origin[2] + spacing[2] * z;
    /// DEBUG ///
    if (curr_z > 1.624) {
      break;
    }

    for (int y = 0; y < ny; y++) {
      int index = (z * ny + y) * nx;
      double curr_label = fine_label_field->GetPointData()->GetScalars()->GetTuple1(index);
      if (curr_label == point_labels[0]
          || curr_label == point_labels[1]
          || curr_label == point_labels[2]
          || curr_label == point_labels[3]) {
        double dist = square(y * spacing[1] + origin[1]);

        if (min_dist < 0.0 || dist < min_dist) {
          min_dist = dist;
        }
        if (max_dist < dist) {
          max_dist = dist;
        }
      }
    }

    if (min_dist < 0.0) {
      continue;
    }

    for (int x = 0; x < nx; x++) {
      for (int y = 0; y < ny; y++) {
        double dist = square(x * spacing[0] + origin[0]) + square(y * spacing[1] + origin[1]);
        if (min_dist <= dist && dist <= max_dist + 0.1) {
          int index = (z * ny + y) * nx + x;

          if (coarse_label_field->GetPointData()->GetScalars()->GetTuple1(index) == original_label) {
            coarse_label_field->GetPointData()->GetScalars()->SetTuple1(index, num_labels);
          }
        }
      }
    }
  }

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetInputData(coarse_label_field);
  writer->SetFileName(kDigOutputFile);
  writer->Write();

  printf("} dig_test\n");
}

// The minimum x is zero, and we want to create a negative x mirror.
void mirror_test() {
  printf("mirror_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(kScalarFieldFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> mesh =
      vtkSmartPointer<vtkStructuredPoints>::New();
  mesh->DeepCopy(reader->GetOutput());

  int dimensions[3];
  mesh->GetDimensions(dimensions);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  double spacing[3], origin[3];
  mesh->GetSpacing(spacing);
  mesh->GetOrigin(origin);

  printf("dimensions: %d, %d, %d\n", dimensions[0], dimensions[1], dimensions[2]);
  printf("spacing: %lf, %lf, %lf\n", spacing[0], spacing[1], spacing[2]);
  printf("origin: %lf, %lf, %lf\n", origin[0], origin[1], origin[2]);

  int new_dimensions[3] = {dimensions[0] * 2 - 1,
                           dimensions[1],
                           dimensions[2]};
  double new_origin[3] = {origin[0] - spacing[0] * (nx - 1),
                          origin[1],
                          origin[2]};

  vtkSmartPointer<vtkDoubleArray> scalar_array =
      vtkSmartPointer<vtkDoubleArray>::New();
  scalar_array->SetNumberOfComponents(1);
  scalar_array->SetNumberOfTuples(
      new_dimensions[0] * new_dimensions[1] * new_dimensions[2]);

  for (int x = 0; x < new_dimensions[0]; x++) {
    for (int y = 0; y < new_dimensions[1]; y++) {
      for (int z = 0; z < new_dimensions[2]; z++) {
        int new_index[] = {x, y, z};
        int new_code = code(new_index, new_dimensions);

        int old_index[] = {x, y, z};
        if (x >= nx - 1) {
          old_index[0] = x - nx + 1;
        } else {
          old_index[0] = nx - 1 - x;
        }

        int old_code = code(old_index, dimensions);
        double scalar = mesh->GetPointData()->GetScalars()->GetTuple1(old_code);
        scalar_array->SetTuple1(new_code, scalar);
      }
    }
  }

  vtkSmartPointer<vtkStructuredPoints> new_mesh =
      vtkSmartPointer<vtkStructuredPoints>::New();
  new_mesh->SetDimensions(new_dimensions);
  new_mesh->SetOrigin(new_origin);
  new_mesh->SetSpacing(spacing);
  new_mesh->GetPointData()->SetScalars(scalar_array);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();
  writer->SetInputData(new_mesh);
  writer->SetFileName(kMirrorOutputFile);
  writer->Write();

  printf("} mirror_test\n");
}

#define BOUNDARY

int main() {
  // threshold_test();
  // waterproof_test();

#ifndef BOUNDARY
  component_test();
  inflate_test();
  deflate_test();
  mask_test();
#else
  // boundary_test();
#endif

  // Specifically for convective cell.
  // crop_test();
  // symmetricalize_test();
  // modify_test();
  // dig_test();
  mirror_test();

  return 0;
}

