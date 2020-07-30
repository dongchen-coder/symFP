#include "sampleNumAnalysis.hpp"

using namespace llvm;

static cl::opt<double>
samplingRateEachLoop("spsrate", cl::Hidden, cl::init(0.0),
                     cl::desc("The sampling rate for each loop"));

namespace sampleNumAnalysis {
    
    char SampleNumberAnalysis::ID = 0;
    static RegisterPass<SampleNumberAnalysis> X("sampleNumAnalysis", "Sample number analysis Pass", false, false);
    
    SampleNumberAnalysis::SampleNumberAnalysis() : FunctionPass(ID) {}
    
    string SampleNumberAnalysis::getBound(Value *bound) {
        
        if (isa<Instruction>(bound)) {
            
            Instruction *inst = cast<Instruction>(bound);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return "(" + getBound(inst->getOperand(0)) + " + " + getBound(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Sub:
                    return "(" + getBound(inst->getOperand(0)) + " - " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Mul:
                    return "(" + getBound(inst->getOperand(0)) + " * " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return "(" + getBound(inst->getOperand(0)) + " / " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Load:
                    return inst->getOperand(0)->getName().str();
                    break;
                case Instruction::Alloca:
                    return inst->getName();
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }
    
    int SampleNumberAnalysis::getConstantStride(Value* expr) {
        
        if (isa<Instruction>(expr)) {
            Instruction *inst = dyn_cast<Instruction>(expr);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    if (isa<Constant>(inst->getOperand(1))) {
                        return dyn_cast<ConstantInt>(inst->getOperand(1))->getValue().getSExtValue();
                    }
                    break;
                default:
                    break;
            }
        }
        
        return 0;
    }

    int SampleNumberAnalysis::valueOfExpression(llvm::Value *expr, std::map<Value *, int> counter) {
        
//        errs() << "In valueOfExpression: ";
//        for (std::map<Value *, int>::iterator it = counter.begin(), eit = counter.end(); it != eit; ++it) {
            
//            errs() << it->first << " " << it->second << " ";
//        }
//        errs() << "\n";
        
        
        
        if (isa<Instruction>(expr)) {
            
            Instruction *inst = cast<Instruction>(expr);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    
//                    inst->getOperand(0)->dump();
                    
//                    errs() << valueOfExpression(inst->getOperand(0), counter) << " " << valueOfExpression(inst->getOperand(1), counter) << "\n";
                    
                    return  valueOfExpression(inst->getOperand(0), counter) + valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Sub:
                    return valueOfExpression(inst->getOperand(0), counter) - valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Mul:
                    return  valueOfExpression(inst->getOperand(0), counter) * valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return valueOfExpression(inst->getOperand(0), counter) / valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Load:
                    return valueOfExpression(inst->getOperand(0), counter);
                    break;
                case Instruction::Alloca:
                    return counter[expr];
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(expr)) {
            return dyn_cast<ConstantInt>(expr)->getValue().getSExtValue();
        }
        
        return 0;
    }
    
    void SampleNumberAnalysis::calculateSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
        if (LoopRefTree->L != NULL) {
            /* calculate sample number */
            
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops_New = loops;
            loops_New.push_back(LoopRefTree);
            
            bool allConstantBounds = true;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops_New.begin(), elit = loops_New.end(); lit != elit; ++lit) {
                for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                    std::string s_tmp = getBound((*(*lit)->LIS->LB)[i].first);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            allConstantBounds = false;
                            break;
                        }
                    }
                    if (allConstantBounds == false) {
                        break;
                    }
                    
                    s_tmp = getBound((*(*lit)->LIS->LB)[i].second);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            allConstantBounds = false;
                            break;
                        }
                    }
                    if (allConstantBounds == false) {
                        break;
                    }
                }
                if (allConstantBounds == false) {
                    break;
                }
            }
            
            
            if (allConstantBounds == true) {
                /* all constant bounds, calculate the number */

//                errs() << "all constants \n";
                
                uint64_t sampling_num = 1;
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops_New.begin(), elit = loops_New.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        int stride = getConstantStride((*(*lit)->LIS->INC)[i]);
                        if (stride < 0) {
                            if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_UGE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SGE) {
                                sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second)) + 1) / (uint64_t) (-stride);
                            } else {
                                sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second))) / (uint64_t) (-stride);
                            }
                        } else {
                            if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_ULE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SLE) {
                                sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first)) + 1) / (uint64_t) (stride);
                            } else {
                                sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first))) / (uint64_t) (stride);
                            }
                        }
                    }
                }

                iter_num += sampling_num;
                
                for (int i = 0; i < loops_New.size(); i++) {
                    sampling_num = sampling_num * samplingRate;
                }
                
                sampleNum[LoopRefTree] = sampling_num;
                
//                errs() << "sample num (calculated) " << sampling_num << "\n";
                
            } else {
                
                std::map<Value*, int> counter;
                std::map<int, Value*> counter_idx;
                
                vector<bool> lbConstant(loops_New.size());
                vector<bool> ubConstant(loops_New.size());
                
                /* init counter idx */
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    if (isa<LoadInst>((*loops_New[i]->LIS->IDV)[0])) {
                        LoadInst* tmp = dyn_cast<LoadInst>((*loops_New[i]->LIS->IDV)[0]);
                        counter_idx[i] = tmp->getOperand(0);
                    } else if (isa<AllocaInst>((*loops_New[i]->LIS->IDV)[0])) {
                        counter_idx[i] = (*loops_New[i]->LIS->IDV)[0];
                    } else {
                        errs() << "Error in init counter\n";
                    }
//                    errs() << counter_idx[i] << "\n";
                }
                
                /* check constant bounds */
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    lbConstant[i] = true;
                    ubConstant[i] = true;
                    
                    std::string s_tmp = getBound((*(loops_New[i])->LIS->LB)[0].first);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            lbConstant[i] = false;
                            break;
                        }
                    }
                    
                    s_tmp = getBound((*(loops_New[i])->LIS->LB)[0].second);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            ubConstant[i] = false;
                            break;
                        }
                    }
                }

                /* init counter */
                errs() << "init counter: ";
                
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    if (i == 0) {
                        counter[counter_idx[i]] = stoi(getBound((*loops_New[i]->LIS->LB)[0].first));
                    } else {
                        if (lbConstant[i] == true) {
                            counter[counter_idx[i]] = stoi(getBound((*loops_New[i]->LIS->LB)[0].first));
                        } else {
                            counter[counter_idx[i]] = valueOfExpression((*loops_New[i]->LIS->LB)[0].first, counter);
                            //(*loops_New[i]->LIS->IDV)[0]->dump();
                        }
                    }
                    
                    errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
                    
                }
                errs() << "\n";
                
                errs() << "Dump stride: ";
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    errs() << getConstantStride((*loops_New[i]->LIS->INC)[0]) << " ";
                }
                errs() << "\n";
                
                /* iterate to get sample number */
                
                uint64_t sampling_num = 0;
                int i_idx = loops_New.size() - 1;
                
//                errs() << "num of loops: " << loops_New.size();
//                LoopRefTree->L->dump();
                
                vector<int> stride(loops_New.size());
                
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    stride[i] = getConstantStride((*loops_New[i]->LIS->INC)[0]);
                }
                
                while(1) {
                    
/*
                    errs() << "counter: ";
                    for (unsigned long i = 0; i < loops_New.size(); i++) {
                        errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
                    }
//                    errs() << sampling_num;
                    errs() << "\n";
*/
 
//                    (*loops_New[i_idx]->LIS->LB)[0].second->dump();
                    
                    int last_level_ub = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter);
                    
//                    errs() << "last level: " << last_level_ub << " - " << counter[(*loops_New[i_idx]->LIS->IDV)[0]] << "\n";
                    if (stride[i_idx] < 0) {
                
                        int last_level_lb = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].first, counter);
                        if (last_level_lb > last_level_ub) {
                            sampling_num += last_level_lb - last_level_ub;
                        }
                        if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_UGE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SGE) && last_level_lb == last_level_ub) {
                            sampling_num += 1;
                        }
                        
                    } else {
                        
                        int last_level_lb = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].first, counter);
                        
                        if (last_level_lb < last_level_ub) {
                            sampling_num += last_level_ub - last_level_lb;
                            
//                            errs() << getBound((*loops_New[i_idx]->LIS->LB)[0].first) << "\n";
//                            errs() << "last level +: " << last_level_ub << " " << last_level_lb << " " << last_level_ub - last_level_lb << "\n";
                            
                        }
                        if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_ULE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SLE) && last_level_lb == last_level_ub) {
                            sampling_num += 1;
                        }
                    }
                    
                    
                    /* here */
//                    errs() << "--------------\n";
                    while(i_idx - 1 >= 0) {
                        i_idx -= 1;
                        
//                        counter[(*loops_New[i_idx]->LIS->IDV)[0]] += stride[i_idx];
                        counter[counter_idx[i_idx]] += stride[i_idx];
                        
                        if (stride[i_idx] < 0) {
                            
                            if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_UGE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SGE)) {
//                                if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] >= valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                if (counter[counter_idx[i_idx]] >= valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                    break;
                                }
                            } else {
//                                if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] > valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                if (counter[counter_idx[i_idx]] > valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                    break;
                                }
                            }
                            
                            
                        } else {
                            
//                            errs() << "Bound: " << getBound((*loops_New[i_idx]->LIS->LB)[0].second) << " ";
//                            errs() << "compare: " << counter[counter_idx[i_idx]] << " " << valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter) << " idx: " << i_idx <<"\n";
                            
                            if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_ULE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SLE)) {
//                                if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] <= valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                if (counter[counter_idx[i_idx]] <= valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                    break;
                                }
                            } else {
//                                if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] < valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                if (counter[counter_idx[i_idx]] < valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                                    break;
                                }
                            }
                        }
                    }
/*
                    errs() << "------------\n";
                    errs() << "counter: ";
                    for (unsigned long i = 0; i < loops_New.size(); i++) {
                        errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
                    }
                    //                    errs() << sampling_num;
                    errs() << "\n";
                    errs() <<"=========";
*/
                    if (i_idx == 0 ) {
                        if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_ULE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SLE)) {
                            
//                            if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] > valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                            if (counter[counter_idx[i_idx]] > valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                                break;
                            
                        } else if (((*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_UGE || (*loops_New[i_idx]->LIS->PREDICATE)[0] == llvm::ICmpInst::ICMP_SGE)) {
                            
//                            if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] < valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                            if (counter[counter_idx[i_idx]] < valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                                break;
                            
                        } else {
                            
//                            if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] == valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                            if (counter[counter_idx[i_idx]] == valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter))
                                break;
                            
                        }
                    }
                    
                    i_idx++;
           
//                    errs() << " " << i_idx << "\n";
                    
                    
                    while(i_idx < loops_New.size()) {
//                        counter[(*loops_New[i_idx]->LIS->IDV)[0]] = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].first, counter);
                        counter[counter_idx[i_idx]] = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].first, counter);
                        i_idx++;
                    }
/*
                    errs() << "counter: ";
                    for (unsigned long i = 0; i < loops_New.size(); i++) {
                        errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
                    }
                    //                    errs() << sampling_num;
                    errs() << "\n\n\n";
*/
                    
                    
                    i_idx = loops_New.size() - 1;
                }
                
//                for (unsigned long i = 0; i < loops_New.size(); i++) {
//                    errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
//                }
//                errs() << "\n";
                
//                errs() << "sample number (executed) " << sampling_num << "\n";
                
//                exit(0);

                for (int i = 0; i < loops_New.size(); i++) {
                    sampling_num = sampling_num * samplingRate;
                }
                
                sampleNum[LoopRefTree] = sampling_num;
                
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops_New = loops;
                if (LoopRefTree->L != NULL) {
                    loops_New.push_back(LoopRefTree);
                }
                calculateSampleNum(*it, loops_New);
            }
        }
        
        return;
    }
    
    void SampleNumberAnalysis::dumpSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LTroot, std::string prefix) {
        
        if (LTroot->L != NULL) {
            if (sampleNum.find(LTroot) != sampleNum.end()) {
                errs() << prefix << "Sample number: " << sampleNum[LTroot] << "\n";
            } else {
                errs() << prefix << "Sample number: N/A \n";
            }
        }
        
        if (LTroot->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
                dumpSampleNum(*it, prefix + "--");
            }
        }
        
        return;
    }
    
    bool SampleNumberAnalysis::runOnFunction(Function &F) {
        
        errs() << " /* Start to analysis the number of samples\n";
        
        samplingRate = samplingRateEachLoop.getValue();
        if (samplingRate == 0.0) {
            samplingRate = 0.01;
        }
        
        /* reading info from previous passes */
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        
        errs() << "calculating:\n";
        iter_num = 0;
        calculateSampleNum(LoopRefTree, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>());
        
        errs() << "Dump tree:\n";
        
        dumpSampleNum(LoopRefTree, "--");
        
        errs() << " End of sample analysis */\n";
        
        return false;
    }
    
    void SampleNumberAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        
        return;
    }
    
}
