#include <femfunc.h>

void gpuCalculateFEM_EbE_vec2(DEVICE_NAME deviceType, dataKeeper &FEMdata, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  ElementsData *elemsData = ElementsData::setElementsData(deviceType, FEMdata);

  elemsData->genMask();
  elemsData->calculateKlocals();

  elemsData->get_Klocals()->getDiagonal(6 * (FEMdata.get_dim() - 1), *elemsData->get_diagK());

//  elemsData->getDiagonalElements(*elemsData->get_Klocals(), *elemsData->get_diagK());
//  elemsData->calculateDiag(*elemsData->get_diagK(), 0.f, 1.f);

  elemsData->initFlocals(0.f, FEMdata.getWaveletParams());

  // (0a)
  elemsData->get_Flocals()->copy(*elemsData->get_r());

  gpuPCG_EbE_vec2(*FEMdata.get_displacements(), deviceType, FEMdata, *elemsData, true, 1e-8f, PRINT_DEBUG_INFO);

  delete elemsData;
}

void gpuPCG_EbE_vec2(Matrix &displacements, DEVICE_NAME devType, dataKeeper &FEMdata, ElementsData &elemsData, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO) {
  CheckRunTime(__func__)
  size_t n_elems  = FEMdata.get_elementsCount();
  size_t DIM = FEMdata.get_dim();
  size_t n_gl_dofs = FEMdata.get_nodesCount() * DIM;
  int grid_size = 6 * (DIM - 1) * n_elems;

//  Matrix *r = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *m = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *z = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *s = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *p = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *u = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *x = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  r->setTo(0.f); m->setTo(0.f); z->setTo(0.f); s->setTo(0.f);
//  p->setTo(0.f); u->setTo(0.f); x->setTo(0.f);

//  elemsData.zeroCGMData();

  Matrix *r = elemsData.get_r(), *m = elemsData.get_m(),
         *z = elemsData.get_z(), *s = elemsData.get_s(),
         *p = elemsData.get_p(), *u = elemsData.get_u(),
         *x = elemsData.get_x();


  // (0a)
//  elemsData.get_Flocals()->copy(*r);
  // (0b)
  elemsData.reductionWithMaskAndTransform(*elemsData.get_diagK(), *m, n_gl_dofs);
  // (0c)
  r->divideElementwise(*m, *z);
  // (0d)
  elemsData.reductionWithMaskAndTransform(*z, *s, n_gl_dofs);
  float gamma0 = r->dotProduct(*s);
  float gamma = gamma0;
  float gamma_new = 0.0f;

  // (0e)
  s->copy(*p);

  if (PRINT_DEBUG_INFO) {
    //std::cout.precision(16);
    std::cout << "gamma0\t\t\t= " << gamma0 << std::endl << std::endl;
  }

  int n_iter = 0;
  do {
    ++n_iter;
    if (PRINT_DEBUG_INFO) {
      std::cout << "Iteration #" << n_iter << std::endl;
    }

    // (1a)
    u->bmm(*p, 1, 6 * (DIM - 1), true,
           *elemsData.get_Klocals(), 6 * (DIM - 1), false, n_elems);

    // (1b)
    float sumElem = p->dotProduct(*u);
    // (1c,d)
    float alpha = gamma / sumElem;
    // (2a)
    x->addWeighted(*p, 1.f, alpha);
    // (2b)
    r->addWeighted(*u, 1.f, -1.f * alpha);
    // (3)
    r->divideElementwise(*m, *z);
    // (4)
    elemsData.reductionWithMaskAndTransform(*z, *s, n_gl_dofs);
    gamma_new = r->dotProduct(*s);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "alpha (gamma / sumElem)\t= " << alpha << std::endl;
      std::cout << "alpha numerator (gamma)\t= " << gamma << std::endl;
      std::cout << "alpha denominator\t= " << sumElem << std::endl;
      // -----------------------------------------------------------------------
    }
    // (5)
    if (gamma_new < eps * gamma0)
      break;

    // (6)
    p->addWeighted(*s, gamma_new / gamma, 1.f);

    if (PRINT_DEBUG_INFO) {
      // Verbose
      // -----------------------------------------------------------------------
      std::cout << "beta\t\t\t= " << gamma_new / gamma << std::endl;
      std::cout << "gamma_new\t\t= " << gamma_new << std::endl;
      std::cout << std::endl;
      // -----------------------------------------------------------------------
    }
    gamma = gamma_new;
    if (!std::isfinite(gamma)) {
      std::cout << "gamma is nan/inf\n";
      exit(-1);
    }
  } while (1);

  if (doAssemblyRes) {
    Matrix *temp = Matrix::setMatrix(elemsData.get_device(), FEMdata.get_nodesCount(), DIM);
    elemsData.reductionWithMask(*x, *temp);
    temp->divideElementwise(*elemsData.get_adjElements(), *temp);
    temp->copy(displacements);
    delete temp;
  } else {
    x->copy(displacements);
  }

//  std::cout << "RESULT:\n";
//  std::cout << FEMdata.get_displacements()->get_numRows() << " " << FEMdata.get_displacements()->get_numCols() << "\n";
//  FEMdata.get_displacements()->Show();


}

void gpuCalculateFEM_DYN2(DEVICE_NAME devType, dataKeeper &dk, bool PRINT_DEBUG_INFO) {
  // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
  CheckRunTime(__func__)
  //assert(beta2 >= beta1 && beta1 >= 0.5f);

  float beta1 = dk.getMechParams().beta1;
  float beta2 = dk.getMechParams().beta2;
  float damping_alpha = dk.getMechParams().damping_alpha;
  float damping_beta = dk.getMechParams().damping_beta;
  float dt = dk.getMechParams().dt;
  float endtime = dk.getMechParams().endtime;

  size_t n_elems = dk.get_elementsCount();
  size_t DIM = dk.get_dim();

  bool isExplicit = beta2 == 0;
  if (isExplicit)
    assert(damping_beta == 0.0f);

  bool doAssemblyRes = false;
  bool isLumped = isExplicit;
  bool isDamping = false;//!(damping_alpha == 0.0f && damping_beta == 0.0f);

  float eps_PCG = 1e-4f;

  float eps_relax;
  bool is_relax;

  Matrix *displ = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  Matrix *vel = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  Matrix *x_dyn1 = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
  displ->setTo(0.f); vel->setTo(0.f); x_dyn1->setTo(0.f);

//  Matrix *x = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  Matrix *r = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
//  x->setTo(0.f); r->setTo(0.f);

//  AssignLoadElement(FEMdata, nodeAdjElem);
//  ApplyLoads_EbE(FEMdata);
//  ApplyConstraints_EbE(FEMdata);

  ElementsData *elemsData = ElementsData::setElementsData(devType, dk);
  elemsData->genMask();
  elemsData->calculateKlocals();
  elemsData->calculateMlocals(isLumped, dk.getMechParams());
//  elemsData->get_Mlocals()->getRow(1)->Show();
//  elemsData->get_Mlocals()->scale(0.5f);
  elemsData->initFlocals(0.f, dk.getWaveletParams());

  Matrix *r = elemsData->get_r();
  Matrix *x_dyn = elemsData->get_x();

//  elemsData->get_Mlocals()->writeToFile("C:/Users/mokin/Desktop/fem_stuff_test/out.txt");
//  elemsData->getDiagonalElements(*elemsData->get_Klocals(), *elemsData->get_diagK());
//  elemsData->get_diagM()->setTo(0.f);
//  elemsData->calculateDiag(*elemsData->get_diagM(), cM, cK, cC, damping_alpha, damping_beta);

  // cM, cK, cC
  float cM = 1.f, cK = 0.f, cC = 0.f;

  Matrix *dampTemp;
  if (isDamping) {
    dampTemp = Matrix::setMatrix(devType, n_elems, 6 * (DIM - 1));
    if (isExplicit) {
      cK = 0.f;
      cC = beta1 * dt;
    } else {
      cK = 0.5f * beta2 * dt * dt;
      cC = beta1 * dt;
    }
  } else {
    if (isExplicit) {
//      elemsData->get_Mlocals()->getDiagonal(6 * (DIM - 1), *elemsData->get_diagM());
    } else {
      cK = 0.5f * beta2 * dt * dt;
      cC = 0.f;

      elemsData->get_Klocals()->getDiagonal(6 * (DIM - 1), *elemsData->get_diagK());
      elemsData->get_Mlocals()->getDiagonal(6 * (DIM - 1), *elemsData->get_diagM());

  //    elemsData->get_Klocals()->scale(cK);

      elemsData->get_diagK()->addWeighted(*elemsData->get_diagM(), cK, cM);
    }
  }


  int endnt;
  is_relax = (endtime < 0.0f);
  if (!is_relax) {
    endnt = static_cast<int>(endtime / dt);
  } else {
    eps_relax = std::fabs(endtime);
    if (PRINT_DEBUG_INFO)
      std::cout << "eps_relax = " << eps_relax << "\n";
  }
/*
  std::string out_path = "C:/Users/mokin/Documents/CAE-Fidesys-5.2/test_segy/test.txt";

  size_t segy_step = 1;
  int axis = 1;
  std::vector<int> receivers = {14,146,2,63,4,55,6,27,81,93};//{130, 78, 146};
  size_t sampleNumber = 0;

  std::fstream out;
  out.open(out_path, std::ios::out);
  out << receivers.size() << "\t";
  for (size_t s = 0; s < receivers.size(); ++s) {
    out << receivers[s] << "\t";
  }
  out << "\n";

  out.close();
*/

  int nt = 1;
  float cnorm_acc, cnorm_vel;
  do {
    float t = nt * dt;
    if (PRINT_DEBUG_INFO) {
      std::cout << "======= Time iteration #" << nt;
      if (!is_relax) std::cout << "/" << endnt << " =======";
      std::cout << " =======";
      std::cout << "\n=========== Time " << t << " ===========\n\n";
    }

//    std::cout << displ->min() << " " << displ->max();

    displ->addWeighted(*vel, 1.f, dt);
    displ->addWeighted(*x_dyn, 1.f, 0.5f * (1.f - beta2) *  dt * dt);
    vel->addWeighted(*x_dyn, 1.f, (1.f - beta1) * dt);
    x_dyn->setTo(0.f);

    r->bmm(*displ, 1, 6 * (DIM - 1), false,
           *elemsData->get_Klocals(), 6 * (DIM - 1), false, n_elems);

    if (isDamping) {
      elemsData->get_Mlocals()->addWeighted(*elemsData->get_Klocals(),
                                            damping_alpha, damping_beta,
                                            *elemsData->get_Clocals());

      dampTemp->bmm(*vel, 1, 6 * (DIM - 1), true,
                    *elemsData->get_Clocals(), 6 * (DIM - 1), false, n_elems);
      r->add(*dampTemp);
    }


    elemsData->calcFlocals(t, dk.getWaveletParams());

    r->addWeighted(*elemsData->get_Flocals(), -1.f, 1.f);

    // TODO: Think how not to use loadVectors! Too dificult!
    /*
    std::unordered_map <int, CPU_Matrix> loadVectors;
    loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                            // instead of assigned. See if-statement in for-loop in the function's body
    GetMapElement2Loadvector(FEMdata, loadVectors, t);
    copyLoads(gpu_data, loadVectors, FEMdata.DIM, n_elems);

    gpuAdd(gpu_data.get_r(), gpu_data.get_loads(), grid_size);
    */

    if (isDamping) {
      if (isExplicit) {
        elemsData->solveDiagSystem(*elemsData->get_diagK(), *r, *x_dyn, doAssemblyRes);
      } else {
//          gpuPCG_EbE_vec_DYN_DAMP(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);
      }
    } else {
      if (isExplicit) {
        elemsData->solveDiagSystem(*elemsData->get_Mlocals(), *r, *x_dyn, doAssemblyRes);
      } else {
        gpuPCG_EbE_vec2(*x_dyn, devType, dk, *elemsData, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);
//        gpuPCG_EbE_vec_DYN(FEMdata, gpu_data, doAssemblyRes, eps_PCG, PRINT_DEBUG_INFO);
      }
    }

    vel->addWeighted(*x_dyn, 1.f, beta1 * dt);
    if (!isExplicit) {
      displ->addWeighted(*x_dyn, 1.f, 0.5f * beta2 * dt * dt);
    }

    ++nt;

    if (is_relax) {
      // TODO: add CNorm !
      cnorm_acc = 0.f;//gpuCNorm(gpu_data.get_x(), grid_size);
      cnorm_vel = 0.f;//gpuCNorm(gpu_data.get_vel(), grid_size);
      if (PRINT_DEBUG_INFO) {
        std::cout << "C-norm acc = " << cnorm_acc << "\nC-norm vel = " << cnorm_vel << std::endl;
        std::cout << std::endl;
      }
      if ((cnorm_vel < eps_relax)) // && (cnorm_acc < eps_relax)
        break;
    } else {
      if (PRINT_DEBUG_INFO) {
        std::cout << std::endl;
      }

      if (nt > endnt) break;
    }

//    if (nt % segy_step == 0) {
//      Matrix *src = Matrix::setMatrix(elemsData->get_device(), dk.get_nodesCount(), DIM);
//      src->setTo(0.f);
//      elemsData->reductionWithMask(*vel, *src);
//      src->divideElementwise(*elemsData->get_adjElements(), *src);

//      writeDisplForSEGY(out_path, *src, receivers, axis, nt * dt);
//      sampleNumber++;

//      delete src;
//    }
  } while (true);

//  std::string path1 = "C:/Users/mokin/Documents/CAE-Fidesys-5.2/SegyTest/fidesys01/fidesys01_Vy.sgy";
//  std::string path2 = "C:/Users/mokin/Documents/CAE-Fidesys-5.2/SegyTest/fidesys01/fidesys01_Vy_TEST.sgy";

//  testSEGY();

//  std::string fidesys_segy = "C:/Users/mokin/Documents/CAE-Fidesys-5.2/SegyTest/fidesys01/fidesys01_Ax.sgy";
//  getSEGYinfo(fidesys_segy);

  std::string SEGYpath = "C:/Users/mokin/Documents/CAE-Fidesys-5.2/test_segy/test.sgy";
//  SegY segy1(0);

//  segy1.addFilenameTXT(out_path, "");
//  segy1.setParameter(SEGYpath, ReceiversType::DISPLACEMENT, 2, *dk.get_nodes(), n_elems, dk.get_nodesCount(), nt, true);

//  convert(SEGYpath, out_path, ReceiversType::VELOCITY, *dk.get_nodes(), 2, dk.get_elementsCount(), dk.get_nodesCount());

//  convertToSEGY(out_path, SEGYpath, receivers, dt, sampleNumber);
//  getSEGYinfo(SEGYpath);


  // Prepare final displacements
  Matrix *temp = Matrix::setMatrix(elemsData->get_device(), dk.get_nodesCount(), DIM);
  elemsData->reductionWithMask(*displ, *temp);
  temp->divideElementwise(*elemsData->get_adjElements(), *temp);
  temp->copy(*dk.get_displacements());
  delete temp;


  delete elemsData;

  delete displ;
  delete vel;
  delete x_dyn1;

//  if (dampTemp)
//    delete dampTemp;
}


void WriteResults(dataKeeper &dk, ResultsDataKeeper &RESdata) {
  if (dk.get_dim() == 2) {
    MakeVTKfile2D(dk.getDataPaths().output_vtk, *dk.get_nodes(), *dk.get_elementsIds(), *dk.get_displacements());
  } else if (dk.get_dim() == 3) {
    MakeVTKfile3D(dk.getDataPaths().output_vtk, *dk.get_nodes(), *dk.get_elementsIds(), *dk.get_displacements());
  }
}

void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO) {

  if (FEMdata.DIM == 2) {
    MakeVTKfile2D(output_vtk, FEMdata.nodes, FEMdata.elements,
                *FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises, RESdata.SmoothStress);
  } else if (FEMdata.DIM == 3) {
    MakeVTKfile3D(output_vtk, FEMdata.nodes, FEMdata.elements,
                *FEMdata.displacements, RESdata.Stress, RESdata.sigma_mises, RESdata.Deformation, RESdata.epsilon_mises);
  }
}

